An example of using **Ant Colony Algorithm** for **task scheduling problem in Mobile Robot Fulfillment Systems**.

------

## Mathematical Models

A mixed integer programming model is presented to describe the problem  as follows:

<img src="Pictures/Parameters.png" style="zoom:67%;" />



<img src="Pictures/model.png" style="zoom: 67%;" />

------

## Heuristic Solutions

- ### The Longest Setup & Handling Time First Rule 

The LSHT rule assigns at t = 0 the R largest tasks to the R robots. After that, whenever a robot is freed, the longest
task among those not yet allocated is put on the robot. This heuristic rule attempts to put the shorter tasks towards the
end of the schedule, where they can be used for workload balancing.   

- ### The Shortest Setup Time First Rule  

Similar to the LSHT rule, the SST rule first assigns the R tasks requiring the shortest setup times to the R robots.
Afterwards, whenever a robot is idle and available, the task among those not yet allocated requiring the smallest setup
time is scheduled on the robot. This rule tries to reduce the total setup time (i.e., the total travel time without loads) of
each robot.  

- ### Ant Colony Algorithm

For the multi-robot scheduling problem under discussion, a complete schedule can be constructed by iteratively dispatching a task unassigned to one of the robots. Thus, an ACO algorithm is proposed here to solve this problem.  The algorithm generally follows the scheme of the any colony system (ACS) . The combined heuristic information for a dispatch and the pheromone updating rules  can be found in the [paper](https://ieeexplore.ieee.org/document/9177514). The C++ code is [here](https://github.com/duohedounai/AntColonyOptimization/blob/main/ACO.cpp).

------

## Numerical Experiments

This paper explicitly studies a multi-robot scheduling problem with the objective of makespan minimization in a
MRFS environment for the first time. Two **heuristic rules (LHST, SST)** and an **ACO algorithm** are developed to solve the problem.   

<img src="Pictures/results.png" style="zoom: 50%;" />

The results show that in most cases all the three proposed heuristics can result in feasible solutions within acceptable time durations, and ACO can generally produce the best schedules.  

------

More details can be found in the paper "**A Task Scheduling Problem in Mobile Robot Fulfillment Systems**". Download the paper at the following website: https://ieeexplore.ieee.org/document/9177514.