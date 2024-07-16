#include "Arduino.h"
#include "MotorControl.h"

// Initialize static variables
volatile int MotorControl::encoderCounts = 0;
volatile unsigned long MotorControl::lastMicrosENCTick = 0;
volatile unsigned long MotorControl::encoderDTmicros = 10000000;
volatile unsigned long MotorControl::encoderDTmicrosFiltered = 10000000;
volatile int MotorControl::rotationDir = 1;

MotorControl::MotorControl(int phasePin, int pwmPin, float supplyVoltage, float gearRatio)
    : _phasePin(phasePin), _pwmPin(pwmPin), _supplyVoltage(supplyVoltage) , _gearRatio(gearRatio) {}

void MotorControl::Init() {
    pinMode(ENCODER_PINA, INPUT);
    pinMode(ENCODER_PINB, INPUT);
    pinMode(_phasePin, OUTPUT);
    pinMode(_pwmPin, OUTPUT);
    digitalWrite(_phasePin, HIGH);
    analogWrite(_pwmPin, 0);

    // Attach interrupt to encoder pins
    attachInterrupt(digitalPinToInterrupt(ENCODER_PINA), MotorControl::encoderA, CHANGE);
    attachInterrupt(digitalPinToInterrupt(ENCODER_PINB), MotorControl::encoderB, CHANGE);
}

// RPM / Rads measured from time between pulses
float MotorControl::ReadRPM() {
    if (encoderDTmicrosFiltered == 0) encoderDTmicrosFiltered = 1;
    unsigned long calculatedEncoderDT = encoderDTmicrosFiltered;
    if ((micros() - lastMicrosENCTick) > MINIMUM_ENCDT) {
        calculatedEncoderDT += (micros() - lastMicrosENCTick);
    }
    float rpmEncoder = rotationDir * ENCDT_RPM / float(calculatedEncoderDT);
    float rpmMotor = rpmEncoder / _gearRatio;
    return rpmMotor;
}

float MotorControl::ReadRads() {
    if (encoderDTmicrosFiltered == 0) encoderDTmicrosFiltered = 1;
    unsigned long calculatedEncoderDT = encoderDTmicrosFiltered;
    if ((micros() - lastMicrosENCTick) > MINIMUM_ENCDT) {
        calculatedEncoderDT += (micros() - lastMicrosENCTick);
    }
    float radsEncoder = rotationDir * ENCDT_RADS / float(calculatedEncoderDT);
    float radsMotor = radsEncoder / _gearRatio;
    return radsMotor;
}

// RPM / Rads measured from number of pulses between function calls
float MotorControl::readRPMPulsesr(){
    static unsigned long lastMillisEnc = 0;
    static unsigned int lastEncoderPos = 0;
    static unsigned int encoderPos = 0;
    unsigned long dtMillisEnc = (millis() - lastMillisEnc);
    if (dtMillisEnc == 0) dtMillisEnc = 1;
    lastMillisEnc = millis(); 
    lastEncoderPos = encoderPos;
    encoderPos = encoderCounts;
    int encoderPulses = encoderPos - lastEncoderPos;
    float rpmEncoder =  ENCPULSES_RPM * (float(encoderPulses) / float(dtMillisEnc));
    float rpmMotor = rpmEncoder / _gearRatio;
    return rpmMotor;
} 

float MotorControl::readRadsPulses(){
    static unsigned long lastMillisEnc = 0;
    static unsigned int lastEncoderPos = 0;
    static unsigned int encoderPos = 0;
    unsigned long dtMillisEnc = (millis() - lastMillisEnc);
    if (dtMillisEnc == 0) dtMillisEnc = 1;
    lastMillisEnc = millis(); 
    lastEncoderPos = encoderPos;
    encoderPos = encoderCounts;
    int encoderPulses = encoderPos - lastEncoderPos;
    float radsEncoder = ENCPULSES_RADS * (float(encoderPulses) / float(dtMillisEnc));
    float radsMotor = radsEncoder / _gearRatio;
    return radsMotor;
}

void MotorControl::motorCommand(float cntrlValue) {
    cntrlValue = cntrlValue / _supplyVoltage;
    if (cntrlValue > 1) cntrlValue = 1;
    if (cntrlValue < -1) cntrlValue = -1;
    if (cntrlValue >= 0) {
        analogWrite(_pwmPin, int(cntrlValue * 255));
        digitalWrite(_phasePin, LOW);
    } else {
        analogWrite(_pwmPin, int(-cntrlValue * 255));
        digitalWrite(_phasePin, HIGH);
    }
}

float MotorControl::StepCMD(float stepPeriod) {
    unsigned long stepPeriodMillis = (unsigned long)(stepPeriod * 1000);
    float stepCmdValue = 0;
    if ((millis() % (stepPeriodMillis * 2)) >= stepPeriodMillis) stepCmdValue = 1;
    return stepCmdValue;
}

float MotorControl::RampCMD(float rampPeriod) {
    unsigned long stepPeriodMillis = (rampPeriod * 1000);
    float rampCmdValue = 0;
    unsigned long rampMillis = millis() % stepPeriodMillis;
    rampCmdValue = float(rampMillis) / (rampPeriod * 1000);
    return rampCmdValue;
}

float MotorControl::SineCMD(float f_Hz) {
    float w_rads = f_Hz * TWO_PI;
    float sineWave = (sin(w_rads * float(millis()) / 1000.0) + 1) / 2;
    return sineWave;
}

void MotorControl::encoderA() {
    #ifndef FILTER_PPR
        EncoderDTMicrosSample();
    #endif
    if (digitalRead(ENCODER_PINA) == HIGH) {
        #ifdef FILTER_PPR
            EncoderDTMicrosSamplePPR();
        #endif
        if (digitalRead(ENCODER_PINB)) {
            rotationDir = 1;
            encoderCounts++;
        } else {
            rotationDir = -1;
            encoderCounts--;
        }
    } else {
        digitalRead(ENCODER_PINB) ? encoderCounts-- : encoderCounts++;
    }
}

void MotorControl::encoderB() {
    #ifndef FILTER_PPR
        EncoderDTMicrosSample();
    #endif
    if (digitalRead(ENCODER_PINB) == HIGH) {
        digitalRead(ENCODER_PINA) ? encoderCounts-- : encoderCounts++;
    } else {
        digitalRead(ENCODER_PINA) ? encoderCounts++ : encoderCounts--;
    }
}

void MotorControl::EncoderDTMicrosSamplePPR() {
    unsigned long microsTick = micros();
    encoderDTmicros = microsTick - lastMicrosENCTick;
    lastMicrosENCTick = microsTick;
    static unsigned long lastENCdt[3] = {10000};
    lastENCdt[0] = lastENCdt[1];
    lastENCdt[1] = lastENCdt[2];
    lastENCdt[2] = encoderDTmicros;
    encoderDTmicrosFiltered = (lastENCdt[0] + lastENCdt[1] + lastENCdt[2]) / 3;
}

void MotorControl::EncoderDTMicrosSample() {
    unsigned long microsTick = micros();
    encoderDTmicros = microsTick - lastMicrosENCTick;
    lastMicrosENCTick = microsTick;
    static unsigned long lastENCdt[ENCODER_FILTER] = {10000};
    unsigned long calcEncDT = encoderDTmicros;
    for (int ii = 0; ii < (ENCODER_FILTER - 1); ii++) {
        lastENCdt[ii] = lastENCdt[ii + 1];
        calcEncDT += lastENCdt[ii];
    }
    lastENCdt[ENCODER_FILTER - 1] = encoderDTmicros;
    encoderDTmicrosFiltered = calcEncDT / ENCODER_FILTER;
}

