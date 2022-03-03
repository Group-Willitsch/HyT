  
#include "max6675.h" //Die MAX6675 Bibliothek
 
int max6675SO = 8; // Serial Output am PIN 8
int max6675CS = 9; // Chip Select am PIN 9
int max6675CLK = 10; // Serial Clock am PIN 10
 
// Initialisierung der MAX6675 Bibliothek mit 
// den Werten der PINs
MAX6675 ktc(max6675CLK, max6675CS, max6675SO); 
 
  
void setup() {
  Serial.begin(9600); // begin serial with baudrate 9600
  delay(500); // take a break in ms
}
 
void loop() {
  // read it
  Serial.print(ktc.readCelsius());
  Serial.println("C"); 
 
 
  // 500ms sleep
  delay(5000);
}
