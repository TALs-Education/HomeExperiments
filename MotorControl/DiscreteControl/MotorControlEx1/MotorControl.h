#ifndef MotorControl_H_
#define MotorControl_H_

#include "Arduino.h"

// define pins
#define ENCODER_PINA 2
#define ENCODER_PINB 3

// define filter mode PPR (3 pulses per revolution - less controller resources), for CPR mode (12 counts per revolution) - uncomment #define FILTER_PPR
//#define FILTER_PPR

// define constants
#define ENCODER_PPR 3         // encoder pulses per revolution
#define ENCODER_CPR 12        // encoder counts per revolution
#define MINIMUM_ENCDT 72000   // define maximum time between pulses to check if the motor has stopped. 1200 --> 60RPM --> 1 sample per control loop CPR Mode 

// define conversion constants 
#define ENCPULSES_RPM 5000     // 12 pulses per revolution -->60/12 = 5 --> 5*1000, since samples are at milliseconds
#define ENCPULSES_RADS 523.59  // ENCPULSES_RPM*2*pi/60
#define RADS_RPM 9.5493       // convert from rad/s to RPM

// define Encoder dt to RPM conversion coefficients
#ifdef FILTER_PPR
    // coefficients for encoder sample at PPR 3 pulses
    #define ENCDT_RPM 20000000 
    #define ENCDT_RADS 2094392 
    #define ENCODER_FILTER 3     // define the number of values to filter for the encoder measurement.
#else
    // coefficients for encoder sample at CPR 12 pulses
    #define ENCDT_RPM 5000000 
    #define ENCDT_RADS 523598 
    #define ENCODER_FILTER 12     // define the number of values to filter for the encoder measurement.
#endif

class MotorControl {
public:
    MotorControl(int phasePin, int pwmPin, float supplyVoltage, float gearRatio);

    // initialize encoder, Motor, attach ISR functions
    void Init();

    // read encoder functions based on averaging time between pulses
    float ReadRPM();
    float ReadRads();
    // read encoder function based on number of pulses over time.
    float readRPMPulsesr();
    float readRadsPulses();
    // motor command update
    void motorCommand(float cntrlValue);

    // motion profiles
    float StepCMD(float stepPeriod);
    float RampCMD(float rampPeriod);
    float SineCMD(float f_Hz);

    // Encoder ISR functions - Interrupt Service Routine
    static void encoderA();
    static void encoderB();

private:
    int _phasePin;
    int _pwmPin;
    float _supplyVoltage;
    float _gearRatio;

    // encoder variables
    static volatile int encoderCounts;
    static volatile unsigned long encoderDTmicros;
    static volatile unsigned long encoderDTmicrosFiltered;
    static volatile unsigned long lastMicrosENCTick;
    static volatile int rotationDir;

    // sample encoder time between pulses
    static void EncoderDTMicrosSample();
    // sample time between pulses simplified for resources, for 3 pulses option ppr mode
    static void EncoderDTMicrosSamplePPR();  
};

#endif
