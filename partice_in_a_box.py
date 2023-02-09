from vpython import *

ball = sphere(pos=vec(-5, 0, 0), radius=0.5, color=color.cyan)
wallR = box(pos=vec(5, 0, 0), size=vec(0.2, 12, 12), color=color.green)
wallL = box(pos=vec(-15, 0, 0), size=vec(0.2, 12, 12), color=color.green)
wallU = box(pos=vec((wallL.pos.x + wallR.pos.x) / 2, 6, 0), size=vec((wallR.pos.x - wallL.pos.x), 0.2, 12),
            color=color.blue)
wallD = box(pos=vec((wallL.pos.x + wallR.pos.x) / 2, -6, 0), size=vec((wallR.pos.x - wallL.pos.x), 0.2, 12),
            color=color.blue)
wallB = box(pos=vec((wallL.pos.x + wallR.pos.x) / 2, 0, -6), size=vec((wallR.pos.x - wallL.pos.x), 12, 0.2),
            color=color.red)

ball.velocity = vec(10, 5, 10)
arrow_scale = 0.1
ball_arrow = arrow(pos=ball.pos, axis=ball.velocity * arrow_scale, color=color.yellow)
deltat = 1 / 30
t = 0

while True:
    # scene.autoscale=False -> Does not autoscale
    rate(1 / deltat)
    if (abs(ball.pos.x - wallR.pos.x) <= ball.radius + wallR.size.x / 2 or
            abs(ball.pos.x - wallL.pos.x) <= ball.radius + wallL.size.x / 2):
        ball.velocity.x = -ball.velocity.x

    if (abs(ball.pos.y - wallU.pos.y) <= ball.radius + wallU.size.y / 2 or
            abs(ball.pos.y - wallD.pos.y) <= ball.radius + wallD.size.y / 2):
        ball.velocity.y = -ball.velocity.y

    if (abs(ball.pos.z - wallB.pos.z) <= ball.radius + wallB.size.z / 2 or
            abs(ball.pos.z - 6) <= ball.radius + wallB.size.z / 2):
        ball.velocity.z = -ball.velocity.z

    ball.pos += ball.velocity * deltat
    ball_arrow.pos = ball.pos
    ball_arrow.axis = ball.velocity * arrow_scale
    t += deltat
