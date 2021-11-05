import RPi.GPIO as GPIO
import time

GPIO.setmode(GPIO.BCM)
GPIO.setwarnings(False)

n = 18

GPIO.setup(n, GPIO.OUT, initial = GPIO.LOW)

for i in range(5):
    GPIO.output(n, GPIO.HIGH)

    time.sleep(0.5)

    GPIO.output(n, GPIO.LOW)

    time.sleep(0.5)