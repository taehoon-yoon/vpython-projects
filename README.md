# vpython-projects
### Toy projects with `vpython` for computational physics
Check out details about vpython https://vpython.org/

## Basic Operation 

- Rotating the camera view: `ctrl + drag left mouse button` or `drag right mouse button`
- Zoom in and out: `scroll wheel` or `alt + drag left mouse button`
- Pan(move up, down, left, right): `shift + drag left mouse button`
- - -
## Mass Spring System
<img src="./IMAGES/spring_network.gif" height="372" width="800">

### Run `spring_network.py`
We implemented mass spring system, where each particles(masses) are connected with each other by spring.
We mainly considered 3 forces to calculate each particle acceleration.

- Spring force
- Gravity
- Damping force

Currently, only `box` objects(where particles can collide. In the above `gif` the green thing, which is comprised of 3 boxes) are supported. But in request I will update other object type such as sphere. (Try sphere! It is easy to implement collision detection algorithm than box object.)

You can also add or remove `box` objects inside the code. Check [line 157](https://github.com/sillsill777/vpython-projects/blob/f5e85b72786f599e6abd5749de5f3304dad52885/spring_network.py#L157) through [line 174](https://github.com/sillsill777/vpython-projects/blob/f5e85b72786f599e6abd5749de5f3304dad52885/spring_network.py#L174) in `spring_network.py`

## Simulating star with SPH(Smoothed Particle Hydrodynamics)

<img src="./IMAGES/SPH_star.gif" height="350" width="350"> &nbsp; &nbsp; &nbsp;
<img src="./IMAGES/SPH_star2.gif" height="350" width="350"> &nbsp; &nbsp; &nbsp;

### Run `SPH.py`

- - -
## References
https://philip-mocz.medium.com/create-your-own-spring-network-simulation-with-python-d9246e4091e5

https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
