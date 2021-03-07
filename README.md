# gl2Body

![Alt Text](https://raw.githubusercontent.com/spastras/gl2Body/main/gl2BodyDemo.gif)

## Description

gl2Body is a simple application I developed in the context of a postgraduate course in Mechanics. It integrates Hamilton's equations of a system consisting of a predifined number of gravitationally interacting bodies and presents the results graphically using OpenGL.

Its original puporse was to simulate the motion of two gravitationally interacting bodies, hence the name "gl2Body". It was later updated so that it works with any number of those.

The main purpose of this application is to provide the user with a visual representation of the evolution of the system using any of four of the most popular integrators.

If you have any suggestions or bug reports please don't hesitate to contact me!

## Integrators

* Euler
* Verlet
* Leapfrog
* Runge-Kutta 4th order

## Dependencies

* libpthread
* libm
* libGL
* libGLU
* libglut

On **Linux**, installing an implementation of OpenGL and GLUT should suffice, e.g. on the latest version of Ubuntu this can be done by executing:
```console
$ sudo apt-get install freeglut3-dev
```

On **Mac OS**, these dependecies can be satisfied by installing Xcode and Apple's Command Line Developer Tools.

## Compilation

You can compile and run gl2Body on Linux or Mac OS by executing:
```console
$ make && ./gl2Body
```

## Known issues

* The current projection method is not ideal for cases in which the center of mass is not at rest
* Initial conditions and options can only be set by editing and recompiling the sources
* The precision of this application is limited due to the use of single precision floating point variables