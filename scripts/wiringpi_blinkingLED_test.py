import wiringpi
import time

wiringpi.wiringPiSetupGpio()

wiringpi.pinMode(18 , 1)


for i in range(5):
    wiringpi.digitalWrite(18, 1) 

    time.sleep(0.1)

    wiringpi.digitalWrite(18, 0)

    #print(wiringpi.digitalRead(18))

    time.sleep(0.1)

wiringpi.softToneCreate(17)

wiringpi.softToneWrite(17, 2250)

time.sleep(4)

wiringpi.digitalWrite(17, 0)
