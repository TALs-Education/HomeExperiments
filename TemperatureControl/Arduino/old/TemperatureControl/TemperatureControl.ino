#include "pwm.h"

// ---------------- Pin configuration ----------------
const int PIN_HEATER_PWM = 9;
const int PIN_LED_PWM    = 10;
const int PIN_FAN_PWM    = 11;

const int PIN_TACH       = 2;

const int PIN_TEMP_AMB   = A0;
const int PIN_TEMP_HEAT  = A1;
const int PIN_TEMP_PLATE = A2;

// ---------------- PWM configuration ----------------
const float PWM_FREQ_HZ = 10000.0f;

// UNO R4 fast PWM objects
PwmOut pwmHeater(PIN_HEATER_PWM);
PwmOut pwmLed(PIN_LED_PWM);
PwmOut pwmFan(PIN_FAN_PWM);

// Duty commands, percent: 0.0 to 100.0
float heaterDuty = 0.0f;
float ledDuty    = 0.0f;
float fanDuty    = 0.0f;

// ---------------- ADC / sensor configuration ----------------
// You said you will use 3.3 V as VREF.
const float ADC_VREF = 3.3f;

// UNO R4 ADC supports higher resolution.
// 14-bit gives values 0..16383.
const int ADC_BITS = 14;
const int ADC_MAX  = (1 << ADC_BITS) - 1;

// MCP9701A transfer function:
// Vout = 0.400 V + 0.0195 V/C * TempC
const float MCP9701_OFFSET_V = 0.400f;
const float MCP9701_SLOPE_V_PER_C = 0.0195f;

// Simple averaging
const int TEMP_SAMPLES = 16;

// ---------------- Fan tach configuration ----------------
// Many 2-wire/3-wire/4-wire PC-style fans output 2 pulses/revolution.
// Change this if your fan datasheet says otherwise.
const float TACH_PULSES_PER_REV = 2.0f;

volatile unsigned long tachPulses = 0;

unsigned long lastTachMillis = 0;
unsigned long lastPlotMillis = 0;

const unsigned long TACH_INTERVAL_MS = 1000;
const unsigned long PLOT_INTERVAL_MS = 250;

float fanRPM = 0.0f;

// ---------------- Serial command buffer ----------------
String commandLine;

// ---------------- Safety limits ----------------
// Adjust these to your experiment.
// Heater will be forced off if heater or plate sensor exceeds this.
const float SAFETY_CUTOFF_C = 95.0f;

// If ADC saturates, temperature reading is no longer trustworthy.
const float SENSOR_SATURATION_C = (ADC_VREF - MCP9701_OFFSET_V) / MCP9701_SLOPE_V_PER_C;


// ---------------- Interrupt service routine ----------------
void tachISR() {
  tachPulses++;
}


// ---------------- Utility functions ----------------
float clampPercent(float x) {
  if (x < 0.0f) return 0.0f;
  if (x > 100.0f) return 100.0f;
  return x;
}


float readVoltageAveraged(int pin) {
  unsigned long sum = 0;

  for (int i = 0; i < TEMP_SAMPLES; i++) {
    sum += analogRead(pin);
    delayMicroseconds(200);
  }

  float adc = (float)sum / (float)TEMP_SAMPLES;
  return adc * ADC_VREF / (float)ADC_MAX;
}


float voltageToMCP9701TempC(float v) {
  return (v - MCP9701_OFFSET_V) / MCP9701_SLOPE_V_PER_C;
}


float readMCP9701TempC(int pin) {
  float v = readVoltageAveraged(pin);
  return voltageToMCP9701TempC(v);
}


void applyOutputs() {
  heaterDuty = clampPercent(heaterDuty);
  ledDuty    = clampPercent(ledDuty);
  fanDuty    = clampPercent(fanDuty);

  pwmHeater.pulse_perc(heaterDuty);
  pwmLed.pulse_perc(ledDuty);
  pwmFan.pulse_perc(fanDuty);
}


// Very small JSON-like parser.
// Accepts commands such as:
// {"heater":40}
// {"fan":80}
// {"led":20}
// {"heater":0,"fan":100,"led":5}
// {"alloff":1}
bool extractJsonNumber(const String& s, const char* key, float& value) {
  String pattern = "\"";
  pattern += key;
  pattern += "\"";

  int keyIndex = s.indexOf(pattern);
  if (keyIndex < 0) return false;

  int colonIndex = s.indexOf(':', keyIndex);
  if (colonIndex < 0) return false;

  int start = colonIndex + 1;

  while (start < s.length() && (s[start] == ' ' || s[start] == '\t')) {
    start++;
  }

  int end = start;
  while (end < s.length()) {
    char c = s[end];
    if ((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.') {
      end++;
    } else {
      break;
    }
  }

  if (end == start) return false;

  value = s.substring(start, end).toFloat();
  return true;
}


void handleCommand(const String& cmd) {
  float value;

  if (extractJsonNumber(cmd, "heater", value)) {
    heaterDuty = clampPercent(value);
  }

  if (extractJsonNumber(cmd, "fan", value)) {
    fanDuty = clampPercent(value);
  }

  if (extractJsonNumber(cmd, "led", value)) {
    ledDuty = clampPercent(value);
  }

  if (extractJsonNumber(cmd, "alloff", value)) {
    if (value != 0.0f) {
      heaterDuty = 0.0f;
      fanDuty = 0.0f;
      ledDuty = 0.0f;
    }
  }

  applyOutputs();
}


void readSerialCommands() {
  while (Serial.available() > 0) {
    char c = Serial.read();

    if (c == '\n' || c == '\r') {
      if (commandLine.length() > 0) {
        handleCommand(commandLine);
        commandLine = "";
      }
    } else {
      commandLine += c;

      if (commandLine.length() > 120) {
        commandLine = "";
      }
    }
  }
}


void updateFanRPM() {
  unsigned long now = millis();

  if (now - lastTachMillis >= TACH_INTERVAL_MS) {
    noInterrupts();
    unsigned long pulses = tachPulses;
    tachPulses = 0;
    interrupts();

    float minutes = (now - lastTachMillis) / 60000.0f;
    fanRPM = ((float)pulses / TACH_PULSES_PER_REV) / minutes;

    lastTachMillis = now;
  }
}


void printPlotterData(float ambC, float heaterC, float plateC) {
  // Arduino Serial Plotter accepts label:value pairs separated by tabs.
  // Keep this as the only regular Serial output while plotting.
  Serial.print("Ambient_C:");
  Serial.print(ambC, 2);

  Serial.print("\tHeater_C:");
  Serial.print(heaterC, 2);

  Serial.print("\tPlate_C:");
  Serial.print(plateC, 2);

  Serial.print("\tFan_RPM:");
  Serial.print(fanRPM, 0);

  Serial.print("\tHeater_PWM:");
  Serial.print(heaterDuty, 1);

  Serial.print("\tFan_PWM:");
  Serial.print(fanDuty, 1);

  Serial.print("\tLED_PWM:");
  Serial.print(ledDuty, 1);

  Serial.println();
}


// ---------------- Arduino setup / loop ----------------
void setup() {
  Serial.begin(115200);

  pinMode(PIN_HEATER_PWM, OUTPUT);
  pinMode(PIN_LED_PWM, OUTPUT);
  pinMode(PIN_FAN_PWM, OUTPUT);

  pinMode(PIN_TACH, INPUT_PULLUP);

  // Use external analog reference on AREF.
  // Connect 3.3 V to AREF before using this.
  analogReference(AR_EXTERNAL);
  analogReadResolution(ADC_BITS);

  // Discard first readings after changing reference/resolution.
  for (int i = 0; i < 10; i++) {
    analogRead(PIN_TEMP_AMB);
    analogRead(PIN_TEMP_HEAT);
    analogRead(PIN_TEMP_PLATE);
    delay(5);
  }

  // Start PWM outputs at 10 kHz and 0%.
  pwmHeater.begin(PWM_FREQ_HZ, 0.0f);
  pwmLed.begin(PWM_FREQ_HZ, 0.0f);
  pwmFan.begin(PWM_FREQ_HZ, 0.0f);

  applyOutputs();

  attachInterrupt(digitalPinToInterrupt(PIN_TACH), tachISR, FALLING);

  lastTachMillis = millis();
  lastPlotMillis = millis();
}


void loop() {
  readSerialCommands();
  updateFanRPM();

  float ambC    = readMCP9701TempC(PIN_TEMP_AMB);
  float heaterC = readMCP9701TempC(PIN_TEMP_HEAT);
  float plateC  = readMCP9701TempC(PIN_TEMP_PLATE);

  // Hard safety cutoff.
  // At 3.3 V reference, MCP9701 saturates around 149 C.
  // Your normal expected range is far below this.
  if (heaterC >= SAFETY_CUTOFF_C || plateC >= SAFETY_CUTOFF_C) {
    heaterDuty = 0.0f;
    ledDuty = 100.0f;   // LED full brightness = fault/status indication
    fanDuty = 100.0f;   // full fan during over-temp
    applyOutputs();
  }

  unsigned long now = millis();

  if (now - lastPlotMillis >= PLOT_INTERVAL_MS) {
    lastPlotMillis = now;
    printPlotterData(ambC, heaterC, plateC);
  }
}