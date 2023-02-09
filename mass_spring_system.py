from vpython import *
import time
import numpy as np

gravitational_acc = 9.8
ball_initial_pos = vec(-0.5, -0.5, 0)
ball_initial_velocity = vec(0, 0, 0)
ball2_initial_pos = vec(0.5, -0.5, 0)
ball2_initial_velocity = vec(0, 0, 0)
ball_exact_initial_pos = vec(0, -0.5, 0)
ball_exact_initial_velocity = vec(0, 0, 0)
scene = canvas(title='mass spring system', width=1700, height=750, center=vector(0, ball_initial_pos.y, 1))
ceiling = box(pos=vec(0, 0, 0), size=vec(2, 0.01, 0.5), color=color.orange)

ball = sphere(pos=ball_initial_pos, radius=0.1, color=color.blue)
ball.mass = 0.5
ball.velocity = ball_initial_velocity
ball.acc_force = vec(0, 0, 0)

ball2 = sphere(pos=ball2_initial_pos, radius=0.1, color=color.red)
ball2.mass = 0.5
ball2.velocity = ball2_initial_velocity
ball2.acc_force = vec(0, 0, 0)

ball_exact = sphere(pos=ball_exact_initial_pos, radius=0.1, color=color.green)
ball_exact.mass = 0.5
ball_exact.velocity = ball_exact_initial_velocity
ball_exact.acc_force = vec(0, 0, 0)

spring = helix(coils=15, radius=0.05, thickness=0.03, pos=vec(ball.pos.x, 0, 0),
               axis=ball.pos - vec(ball.pos.x, 0, 0))
spring.k = 10
spring.L0 = 0.5

spring2 = helix(coils=15, radius=0.05, thickness=0.03, pos=vec(ball2.pos.x, 0, 0),
                axis=ball2.pos - vec(ball2.pos.x, 0, 0))
spring2.k = 10
spring2.L0 = 0.5

spring_exact = helix(coils=15, radius=0.05, thickness=0.03, pos=vec(ball_exact.pos.x, 0, 0),
                     axis=ball_exact.pos - vec(ball_exact.pos.x, 0, 0))
spring_exact.k = 10
spring_exact.L0 = 0.5
'''
m*ddot(x)=-k(x+L0)-mg
x=A*cos(sqrt(k/m)t)-mg/k-L0, omega=sqrt(k/m)t
A=x+mg/k+L0 -> initial condition
'''
A = ball_exact.pos.y + ball_exact.mass * gravitational_acc / spring_exact.k + spring_exact.L0
omega = np.sqrt(spring_exact.k / ball_exact.mass)

label1 = label()
label2 = label()

def hook_force(spring):
    return -spring.k * (spring.length - spring.L0) * spring.axis.norm()


def gravitational_force(object):
    return vec(0, -1, 0) * gravitational_acc * object.mass


dt = 1 / 1000
t = 0
flag = True
while True:
    scene.autoscale = False
    rate(1 / dt)
    t += dt
    if flag is True:
        time.sleep(2)
        flag = False

    ball.acc_force = gravitational_force(ball) + hook_force(spring)
    ball.pos += ball.velocity * dt
    ball.velocity += (ball.acc_force / ball.mass) * dt
    spring.axis = ball.pos - spring.pos

    ball2.acc_force = gravitational_force(ball2) + hook_force(spring2)
    ball2.velocity += (ball2.acc_force / ball2.mass) * dt
    ball2.pos += ball2.velocity * dt
    spring2.axis = ball2.pos - spring2.pos

    ball_exact.pos.y = A * np.cos(omega * t) - ball_exact.mass * gravitational_acc / spring_exact.k - spring_exact.L0
    spring_exact.axis = ball_exact.pos - spring_exact.pos

    label1.pos = ceiling.pos + vector(-0.5, 0.2, 0)
    label1.text = ('Diff : %.5f' % abs(ball.pos.y-ball_exact.pos.y))
    label2.pos = ceiling.pos + vector(0.5, 0.2, 0)
    label2.text = ('Diff : %.5f' % abs(ball2.pos.y-ball_exact.pos.y))