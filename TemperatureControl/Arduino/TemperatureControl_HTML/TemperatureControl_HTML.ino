#include <Arduino.h>
#include "pwm.h"

// ---------------- Pin configuration ----------------
const int PIN_HEATER1_PWM = 9;
const int PIN_HEATER2_PWM = 10;  // optional second heater element
const int PIN_FAN_PWM     = 11;

const int PIN_TACH        = 2;

const int PIN_TEMP_AMB    = A0;

// Beam sensors, spaced along the heating beam from heater1 (D9) to the
// optional heater2 (D10). Beam1 is closest to heater1, Beam5 closest to heater2.
const int NUM_BEAM_SENSORS = 5;
const int PIN_TEMP_BEAM[NUM_BEAM_SENSORS] = { A1, A2, A3, A4, A5 };

// ---------------- PWM configuration ----------------
const float PWM_FREQ_HZ = 10000.0f;

PwmOut pwmHeater1(PIN_HEATER1_PWM);
PwmOut pwmHeater2(PIN_HEATER2_PWM);
PwmOut pwmFan(PIN_FAN_PWM);

float heater1Duty = 0.0f;
float heater2Duty = 0.0f;
float fanDuty      = 0.0f;

// ---------------- ADC / MCP9701 configuration ----------------
// 3.3 V connected to AREF
const float ADC_VREF = 3.3f;

const int ADC_BITS = 14;
const int ADC_MAX  = (1 << ADC_BITS) - 1;

// MCP9701A:
// Vout = 0.400 V + 0.0195 V/C * TempC
const float MCP9701_OFFSET_V = 0.400f;
const float MCP9701_SLOPE_V_PER_C = 0.0195f;

const int TEMP_SAMPLES = 16;

// ---------------- Fan tach configuration ----------------
// Most PC-style fans: 2 pulses per revolution
const float TACH_PULSES_PER_REV = 2.0f;

volatile unsigned long tachPulses = 0;
unsigned long lastTachMillis = 0;
float fanRPM = 0.0f;

// ---------------- Timing ----------------
unsigned long lastTelemetryMillis = 0;
const unsigned long TELEMETRY_INTERVAL_MS = 250;

// ---------------- Safety ----------------
// Both heaters are forced off if any beam sensor exceeds this.
const float SAFETY_CUTOFF_C = 125.0f;
bool faultState = false;

// ---------------- Serial command buffer ----------------
String commandLine;


// ---------------- Interrupt ----------------
void tachISR() {
  tachPulses++;
}


// ---------------- Helpers ----------------
float clampPercent(float x) {
  if (x < 0.0f) return 0.0f;
  if (x > 100.0f) return 100.0f;
  return x;
}


void applyOutputs() {
  heater1Duty = clampPercent(heater1Duty);
  heater2Duty = clampPercent(heater2Duty);
  fanDuty     = clampPercent(fanDuty);

  pwmHeater1.pulse_perc(heater1Duty);
  pwmHeater2.pulse_perc(heater2Duty);
  pwmFan.pulse_perc(fanDuty);
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


// Extracts numeric values from simple JSON commands.
// Examples:
// {"heater":40}
// {"heater2":40}
// {"fan":80}
// {"heater":40,"heater2":0,"fan":80}
// {"alloff":1}
// {"clearfault":1}
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

  if (extractJsonNumber(cmd, "clearfault", value)) {
    if (value != 0.0f) {
      faultState = false;
    }
  }

  if (extractJsonNumber(cmd, "alloff", value)) {
    if (value != 0.0f) {
      heater1Duty = 0.0f;
      heater2Duty = 0.0f;
      fanDuty = 0.0f;
      faultState = false;
      applyOutputs();
      return;
    }
  }

  if (extractJsonNumber(cmd, "fan", value)) {
    fanDuty = clampPercent(value);
  }

  // Do not allow heater restart while fault is active.
  if (extractJsonNumber(cmd, "heater", value)) {
    if (!faultState) {
      heater1Duty = clampPercent(value);
    }
  }

  if (extractJsonNumber(cmd, "heater2", value)) {
    if (!faultState) {
      heater2Duty = clampPercent(value);
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

      if (commandLine.length() > 160) {
        commandLine = "";
      }
    }
  }
}


void updateFanRPM() {
  unsigned long now = millis();

  if (now - lastTachMillis >= 1000) {
    noInterrupts();
    unsigned long pulses = tachPulses;
    tachPulses = 0;
    interrupts();

    float minutes = (now - lastTachMillis) / 60000.0f;
    fanRPM = ((float)pulses / TACH_PULSES_PER_REV) / minutes;

    lastTachMillis = now;
  }
}


void checkSafety(float maxBeamC) {
  if (maxBeamC >= SAFETY_CUTOFF_C) {
    faultState = true;

    heater1Duty = 0.0f;
    heater2Duty = 0.0f;
    fanDuty = 100.0f;

    applyOutputs();
  }
}


void printTelemetryJson(float ambC, const float beamC[NUM_BEAM_SENSORS]) {
  Serial.print("{");

  Serial.print("\"t_ms\":");
  Serial.print(millis());

  Serial.print(",\"ambient_c\":");
  Serial.print(ambC, 2);

  for (int i = 0; i < NUM_BEAM_SENSORS; i++) {
    Serial.print(",\"beam");
    Serial.print(i + 1);
    Serial.print("_c\":");
    Serial.print(beamC[i], 2);
  }

  Serial.print(",\"fan_rpm\":");
  Serial.print(fanRPM, 0);

  Serial.print(",\"heater1_pwm\":");
  Serial.print(heater1Duty, 1);

  Serial.print(",\"heater2_pwm\":");
  Serial.print(heater2Duty, 1);

  Serial.print(",\"fan_pwm\":");
  Serial.print(fanDuty, 1);

  Serial.print(",\"fault\":");
  Serial.print(faultState ? 1 : 0);

  Serial.print("}");
  Serial.println();
}


// ---------------- Setup / Loop ----------------
void setup() {
  Serial.begin(115200);

  pinMode(PIN_HEATER1_PWM, OUTPUT);
  pinMode(PIN_HEATER2_PWM, OUTPUT);
  pinMode(PIN_FAN_PWM, OUTPUT);

  pinMode(PIN_TACH, INPUT_PULLUP);

  // External reference: connect 3.3 V to AREF.
  analogReference(AR_EXTERNAL);
  analogReadResolution(ADC_BITS);

  // Discard initial ADC readings after reference setup.
  for (int i = 0; i < 10; i++) {
    analogRead(PIN_TEMP_AMB);
    for (int j = 0; j < NUM_BEAM_SENSORS; j++) {
      analogRead(PIN_TEMP_BEAM[j]);
    }
    delay(5);
  }

  // Start 10 kHz PWM outputs at 0%.
  pwmHeater1.begin(PWM_FREQ_HZ, 0.0f);
  pwmHeater2.begin(PWM_FREQ_HZ, 0.0f);
  pwmFan.begin(PWM_FREQ_HZ, 0.0f);

  applyOutputs();

  attachInterrupt(digitalPinToInterrupt(PIN_TACH), tachISR, FALLING);

  lastTachMillis = millis();
  lastTelemetryMillis = millis();
}


void loop() {
  readSerialCommands();
  updateFanRPM();

  float ambC = readMCP9701TempC(PIN_TEMP_AMB);

  float beamC[NUM_BEAM_SENSORS];
  float maxBeamC = -1000.0f;
  for (int i = 0; i < NUM_BEAM_SENSORS; i++) {
    beamC[i] = readMCP9701TempC(PIN_TEMP_BEAM[i]);
    if (beamC[i] > maxBeamC) maxBeamC = beamC[i];
  }

  checkSafety(maxBeamC);

  unsigned long now = millis();

  if (now - lastTelemetryMillis >= TELEMETRY_INTERVAL_MS) {
    lastTelemetryMillis = now;
    printTelemetryJson(ambC, beamC);
  }
}
