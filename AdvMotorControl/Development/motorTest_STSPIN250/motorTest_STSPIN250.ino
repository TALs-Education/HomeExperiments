#include <Arduino.h>

// Pin assignments
int phaseControl = 6; // D6 for phase control
int pwmControl = 5;  // D5 for PWM control
int enablePin = 7;   // D7 for enable

int encoderA = 2;    // D2 for encoder
int encoderB = 3;    // D3 for encoder

volatile long encoderPosition = 0; // Position from encoder
long previousPosition = 0;         // Previous position for calculating velocity
unsigned long previousTime = 0;    // Previous time for calculating velocity
float velocity = 0;                // Calculated velocity


// Function prototypes
void encoderAChange();
void encoderBChange();

void setup() {
  // Set motor control pins as outputs
  pinMode(phaseControl, OUTPUT);
  pinMode(pwmControl, OUTPUT);
  pinMode(enablePin, OUTPUT);

  // Set encoder pins as inputs
  pinMode(encoderA, INPUT);
  pinMode(encoderB, INPUT);

    // Attach interrupts for encoder
  attachInterrupt(digitalPinToInterrupt(encoderA), encoderAChange, CHANGE);
  attachInterrupt(digitalPinToInterrupt(encoderB), encoderBChange, CHANGE);

    // Set motor direction and speed
    digitalWrite(enablePin, HIGH); // Enable the motor driver
  // Begin serial communication
  Serial.begin(9600);
}

void loop() {
  if (Serial.available() > 0) {
    int speed = Serial.parseInt(); // Read the speed value from serial input
    if (abs(speed)>0){


    if (speed >= 0) {
      digitalWrite(phaseControl, HIGH); // Set motor direction
    } else {
      digitalWrite(phaseControl, LOW); // Set motor direction
      speed = -speed; // Make speed positive for PWM
    }

    analogWrite(pwmControl, map(speed, 0, 1023, 0, 255)); // Set motor speed
    }
  }
    // Calculate velocity as delta position / delta time
  unsigned long currentTime = millis();
  if (currentTime - previousTime >= 100) { // Update velocity every 100 milliseconds
    long currentPosition = encoderPosition;
    velocity = (float)(currentPosition - previousPosition) / (currentTime - previousTime);
    previousPosition = currentPosition;
    previousTime = currentTime;

    // Output velocity for debugging
      // Optional: Output encoder position for debugging
    Serial.print("Encoder Position: ");
    Serial.println(encoderPosition);
    Serial.print("Velocity: ");
    Serial.println(velocity);
  }
}

// Interrupt service routines for the encoder
void encoderAChange() {
  if (digitalRead(encoderA) == digitalRead(encoderB)) {
    encoderPosition++;
  } else {
    encoderPosition--;
  }
}

void encoderBChange() {
  if (digitalRead(encoderA) != digitalRead(encoderB)) {
    encoderPosition++;
  } else {
    encoderPosition--;
  }
}
