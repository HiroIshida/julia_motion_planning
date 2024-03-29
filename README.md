# UPDATE 2023/09/24
This code has some critical bugs in core algorithm and collision checking. Also the code is quite dirty (I had only a year of programming experience at that time). Please refer the newer implementation https://github.com/HiroIshida/diplan-cpp

# Fast marching tree
This repository currently includes fast marching tree (FMT)[1] implementation in Julia. Specifically, The implementation is for double integrator system [2]. I like FMT algorithm for its simpleness and speed compared to rather complex RRT* algorithm. And, this is kind of my missionary activity for FMT. The GIF below was made by concatenating data generated by this `example.jl` in this repository.

![fig1](https://raw.githubusercontent.com/HiroIshida/julia_motion_planning/master/fig/animation.gif)

I'm a newbie in Julia, so I welcome any suggestion about everything. But please note that even a code written by a newbie is that fast.   

[1] Janson, Lucas, et al. "Fast marching tree: A fast marching sampling-based method for optimal motion planning in many dimensions." The International journal of robotics research 34.7 (2015): 883-921.
[2] Schmerling, Edward, Lucas Janson, and Marco Pavone. "Optimal sampling-based motion planning under differential constraints: the driftless case." Robotics and Automation (ICRA), 2015 IEEE International Conference on. IEEE, 2015.
