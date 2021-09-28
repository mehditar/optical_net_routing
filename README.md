# optical_net_routing
This code simulates an optical network where the nodes are connected via optical fiber cables and each fiber cable has a spectrum distribution as in elastic optical network (EON) architecture.  Also, once the demands for transferring data from nodes to nodes are fed to the network, for each demand, sequentially the code conducts a multi-layer routing, modulation technique and spectrum search, and if proper resources are found the code allocates these resources to the demand until the demand transfers its complete data.  In other words, the code keeps track of the arrival time of the demands, searches for and allocates resources to them and once demands completely served, it releases the associated resources to be used for other demands.  For the routing part, this code includes one heuristic algorithm and one integer linear programming algorithm. For any given run of the code the user can choose to use the heuristic or ILP based algorithm.  The details of the goal of the algorithms, their structures and associated results can be found in my papers at https://ieeexplore.ieee.org/abstract/document/9006849/ and https://ieeexplore.ieee.org/abstract/document/9322185 .


Dependencies:

Python, Pandas, Numpy, Gurobi (for Integer linear programming).

The code has 6 files: 

main.py 

This code imports the main files including the network architecture (the nodes, connections, link lengths, etc.,), the demands (including the arrival times, the bandwidth request, etc.) and also is responsible for tracking the nodes and network resources throughout the network operation (which is implemented through iteration). In this file, another function, demand_updates.py takes care of the resource searches and the resource allocation for each demand

demand_update.py

for any given demand at a time, this function chooses either the heuristic algorithm or the ILP  and searches for resources and based on the search results either allocates the resources, updates the resources or puts the demands in line(the latter,  in case the resource are not found). 

heuristic.py

Tis code contains the heuristic algorithm. This algorithm is responsible for choosing one route, one modulation format and the range of spectrum that the signal will occupy throughout the route. It should be mentioned that for each pair, the code has already access to 50 shortest paths between the pair and it is responsible just for selecting the best based on the network state at that moment. The 50 shortest paths for each pair of the network are separately calculated and stored in a file

bandwidth.py

This is part of the heuristic.py function. For each individual path at any given time, this function searches for the best spectrum range that maximizes the use and efficiency of the spectrum in the network.

integer_linear_programming.py

This codes searches for best route, modulation and spectrum among all possible paths between each pair (not just the 50 shortest paths). 

trans_store_clean.py

This file contains several functions that are responsible for  importing files, storing and cleaning the variables, and transforming the formats of the variables. 

Supplemental files:

1.	a-PATHS.dat, a-INFO.dat, a-PathOrder.dat, a-LinkName.dat: these files contain the first 50 paths between each pair of source and destinations in the network. I used shortest path algorithm (SPA) on MATLAB  to obtain these files.
3.	csvND.dat, csvNL.dat: these file contain the network architecture including the size of the network, the nodes and their connections and the length in kilometers of each link. The specific network, here is US backbone network with 28 nodes and 90 fiber links.  
4.	QN.dat: this file contains all demands and their associated information including the arrival time of the demand, the size of the total data, source and destination nodes for the demand, and transfer completion due time.
5.	Par.dat: contains the network configuration values









