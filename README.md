# CG Engine

This is my implementation of a graphics engine for Stuyvesant High School's MKS66 Computer Graphics course. Everything is written from scratch in C++. 

## Features

- Custom vector and matrix math
- Stack-based structure for relative coordinate systems
- Polygon drawing and rendering for boxes, spheres, toruses
- Phong reflection model for lighting
- Bezier and hermite curves

# Features to implement
- Mesh support for .obj files
- Animation
- Custom/multiple lighting

## Usage
Insert commands into `script` per line

Valid commands:
- `sphere [x] [y] [z] [r]` | render sphere centered at (x,y,z) with radius r
- `torus [x] [y] [z] [r] [R]` | render torus centered at (x,y,z) with thickness r and distance R from center
- `box [x0] [y0] [z0] [h] [w] [d]` | render box with one corner (x0,y0,z0) and height h, width w, depth d
- `push` | add new coordinate system
- `pop` | leave current coord system and go to previous coordinate system
- `clear` | clears screen
- `display` | displays using ImageMagick
- `save [filename]` | saves image as png
- `rotate [x/y/z] [deg]` | rotates coordinate system by deg degrees in either x, y, or z axis
- `scale [a] [b] [c]` | scales each axis by a, b, c
- `move [dx] [dy] [dz]` | translates coordinate system by dx, dy, dz in x,y,z directions

To run the script, run `make` in terminal

Note for WSL, imagemagick might be janky, I used XMing (make sure to run `export DISPLAY=localhost:0` before each session)
