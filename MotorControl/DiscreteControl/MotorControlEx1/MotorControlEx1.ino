/* Motor Control 16/07/2024

   Arduino Uno
   Home Experiment Shield
   
*/
#include "MotorControl.h"

// define constants
#define DT_CONTROL_MILLIS 10  // Control DT suggested range 5-50 mSec --> 200-20 Hz
#define GEAR_RATIO 30         // Motor gear ratio
#define TIMEOUT 10            // Experiment length in seconds
#define SUPPLY_VOLTAGE 5.0    // 5 volts

// define pins
#define PHASE_PIN 4 // when jumper resistor is placed 4-->6 then switch to pin 6
#define PWM_PIN 5

// timeout flag
boolean timeOutFlag = 0;

// Control Variables
float motorRADSdt = 0;
float motorRADSPulses = 0;

// control loop function
void MotorControlLoop(void);
// Control function
float pidControl(float desiredCMD);

// Create an instance of the MotorControl class
MotorControl motorController(PHASE_PIN, PWM_PIN, SUPPLY_VOLTAGE,GEAR_RATIO);

// Init setup
void setup() {
  // initialize encoder and motor control
  motorController.Init();
  pinMode(LED_BUILTIN, OUTPUT);
  
  // init serial
  Serial.begin (115200);
}

void loop() {
  if (!timeOutFlag){
    static unsigned long millisTick = millis();
    if (millisTick >= TIMEOUT*1000) timeOutFlag = 1;
    
    // control loop software implementation
    if ((millis() - millisTick) >= DT_CONTROL_MILLIS){
      millisTick = millis();
      unsigned long microsTick = micros(); // simple loop CPU consumption measurement
      digitalWrite(LED_BUILTIN, HIGH);

      // update Motor RPM variables
      motorRADSdt = motorController.ReadRads(); // rad/s measured from dt between pulses with averaging
      motorRADSPulses = motorController.readRadsPulses(); // rad/s measure from number of pulses over time
      //////////////////////////////
      //////// Student Code ////////
      //////////////////////////////

      // motor control loop function call
      MotorControlLoop();

      //////////////////////////////
      //// End of student Code /////
      //////////////////////////////
      digitalWrite(LED_BUILTIN, LOW);
    }
  } else {
      // turn off motor
      motorController.motorCommand(0);
  }
}

//////////////////////////////
////// Student functions /////
//////////////////////////////

// control loop function
void MotorControlLoop(void){
  // Command profile
  float desiredCMD = 30 * motorController.SineCMD(0.5);  // Sine command, input Hz
  //float desiredCMD = 30 * motorController.RampCMD(10);  // Sine command, input Hz
  
  // motor control function
  float motorCMD = pidControl(desiredCMD);
  
  // update motor, values range of [-1..1]
  motorController.motorCommand(motorCMD);

  // send values to serial plotter
  Serial.print("Desired_[Rad/s]:");
  Serial.print(desiredCMD);
  Serial.print(" , ");
  Serial.print("Error[Rad/s]:");
  Serial.print(desiredCMD - motorRADSdt);
  Serial.print(" , ");
  Serial.print("Motor[Rad/s]_dt:");
  Serial.print(motorRADSdt); // motorRPMPulses
  Serial.print(" , ");
  Serial.print("Motor[Rad/s]_pulses:");
  Serial.println(motorRADSPulses); // motorRPMPulses
} 

// control loop function
float pidControl(float desiredCMD){
  // initialize variables
  static float dt = DT_CONTROL_MILLIS / 1000.0;
  static float ciEr = 0; // integral error
  static float cdEr = 0; // differential error
  static float cEr  = 0; // current error
  static float lEr  = 0; // last error
  
  // control coefficients
  float kp = 0.25;
  float ki = 2.5;
  float kd = 0.001;
  
  // update error
  lEr = cEr;
  cEr = desiredCMD - motorRADSdt;
  ciEr = ciEr + cEr * dt;
  cdEr = (cEr - lEr) / dt;
  
  // clip integrator error
  if (ciEr * ki > 2.5) ciEr = 2.5 / ki;
  if (ciEr * ki < -2.5) ciEr = -2.5 / ki;

  // update control command
  float motorCMD = kp * cEr + ki * ciEr + kd * cdEr;
  return motorCMD;
}

//////////////////////////////
// End of Student functions //
//////////////////////////////
