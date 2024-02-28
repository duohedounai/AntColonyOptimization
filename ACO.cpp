#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<ctime>

/********************Parameter setting*******************/
#define ROBOT_NUM 15
#define TASK_NUM 60
#define POSITION_NUM TASK_NUM  //	 The number of tasks that can be accepted by each AGV
#define T_MAX 50  // The number of iterations
#define ANT_NUM 10
#define NODE_NUM TASK_NUM  // The number of nodes that each ant needs to walk
#define Q0 0.7  // The probability of selecting a binary by the strength of the comprehensive heuristic information
#define P0 0.001 // Default initial pheromone concentration
#define A 1  //importance parameter of tau_{ij}
#define B 3  //importance parameter of s
#define C 1  //importance parameter of u
#define R_L 0.2  // Proportion of local volatilization of pheromones
#define R_G 0.2  // Proportion of global volatilization of pheromones

#define BIG_M 999

using namespace std;

typedef struct ant_node
{
	int combination_num;  // The index of the binary assigned this time
	int task;  // The assignment index
	int AGV;  // The index of the AGV on which the task was performed

}ant_node;

typedef struct avaliable_node
{
	vector<int>combination_num;
	vector<int>task;
	vector<int>AGV;

}availiable_node;

#define random_int(a,b) (rand()%(b-a+1)+a)  // Generate random integers between [a,b].
#define random_double(a,b) (rand()/double(RAND_MAX))  // Generates random floating-point numbers between (0,1).

/************************Function declarations************************/
void ant_path_selection(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM], int current_ant_num, int current_node_num, ant_node(&candidate_set)[ANT_NUM][TASK_NUM * ROBOT_NUM], double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM], double(&finishing_time)[ANT_NUM][ROBOT_NUM], double(&M)[ANT_NUM][ROBOT_NUM], double D[ROBOT_NUM][TASK_NUM], double d[TASK_NUM][TASK_NUM], double h[TASK_NUM], vector<vector<vector<int>>>(&task_of_AGV), vector<int>(&local_update_set), vector<int>(&ant_num_for_local_update_set), int(&last_node_find_flag)[ANT_NUM]);
void local_pheromone_update(double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM], vector<int>(&local_update_set), vector<int>former_local_update_set, int(&last_node_find_flag)[ANT_NUM], vector<int>(&ant_num_for_former_local_update_set), vector<int>(&ant_num_for_local_update_set));
void global_pheromone_update(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM],ant_node(&final_solution)[NODE_NUM],double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM], vector<int>(&local_update_data_all_ant), vector<int>(&least_Cmax_combination_num), double(&Cmax), double(&finishing_time)[ANT_NUM][ROBOT_NUM], vector<vector<vector<int>>>task_of_AGV, int iteration_num, int(&iteration_num_of_Cmax), double(&final_finishing_time)[ROBOT_NUM], vector<int>(&ant_num_of_local_update_data_all_ant));
double Cmax_change_statistics(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM], double(&finishing_time)[ANT_NUM][ROBOT_NUM], int current_ant_num, int current_node_num, int current_AGV, int current_task, double D[ROBOT_NUM][TASK_NUM], double d[TASK_NUM][TASK_NUM], double h[TASK_NUM]);
double AGV_utilization_change_statistics(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM], int current_ant_num, int current_node_num, int current_AGV, int current_task, double D[ROBOT_NUM][TASK_NUM], double d[TASK_NUM][TASK_NUM], double h[TASK_NUM]);

/*************************Main function*************************/
int main()
{
	cout << "ACO starts running！" << endl;
	clock_t start, finish;
	start = clock();
	static double d[TASK_NUM][TASK_NUM] =
	{
		0,5,13,15,8,13,7,12,4,9,14,15,11,4,3,16,17,8,11,17,4,3,11,8,9,4,9,10,7,4,3,3,3,16,11,9,8,2,9,11,3,5,2,18,12,7,7,18,8,12,4,11,1,12,6,11,6,9,10,10,
5,0,8,10,3,8,2,17,5,4,9,10,16,9,6,11,12,11,6,12,7,8,6,3,4,9,4,5,2,9,8,8,8,11,6,10,9,3,14,12,8,10,5,13,7,12,2,13,3,17,7,14,6,7,11,6,1,4,5,5,
13,8,0,16,5,6,6,25,9,12,17,10,24,17,10,11,16,19,14,10,15,10,10,7,4,13,6,5,8,17,16,16,12,17,8,18,17,11,22,20,16,14,13,13,7,20,6,15,5,25,9,22,14,1,13,2,9,4,3,11,
15,10,16,0,11,10,10,9,15,6,1,6,8,17,16,5,2,7,4,6,11,18,6,9,12,19,10,11,8,11,18,16,18,1,8,6,7,13,6,4,12,20,13,3,9,10,12,3,11,9,17,6,16,15,21,14,9,12,13,5,
8,3,5,11,0,5,1,20,8,7,12,7,19,12,9,8,11,14,9,9,10,11,5,2,1,12,1,2,3,12,11,11,11,12,3,13,12,6,17,15,11,13,8,10,4,15,5,10,4,20,10,17,9,4,14,3,4,3,2,6,
13,8,6,10,5,0,6,19,13,6,11,4,18,15,14,5,10,13,8,4,9,16,4,5,4,17,4,3,6,11,16,14,16,11,2,12,11,11,16,14,10,18,11,7,1,14,10,9,9,19,15,16,14,7,19,8,7,8,7,5,
7,2,6,10,1,6,0,19,7,6,11,8,18,11,8,9,10,13,8,10,9,10,4,1,2,11,2,3,2,11,10,10,10,11,4,12,11,5,16,14,10,12,7,11,5,14,4,11,3,19,9,16,8,5,13,4,3,2,3,5,
12,17,25,9,20,19,19,0,16,13,8,15,3,8,15,14,11,6,11,15,10,15,15,18,21,12,19,20,17,8,9,9,13,10,17,7,8,14,3,5,9,11,12,12,18,5,19,12,20,2,16,5,11,24,12,23,16,21,22,14,
4,5,9,15,8,13,7,16,0,9,14,15,15,8,1,16,17,10,11,17,6,3,11,8,9,4,9,10,7,8,7,7,3,16,11,9,8,2,13,11,7,5,4,18,12,11,3,18,4,16,2,13,5,8,6,7,6,5,6,10,
9,4,12,6,7,6,6,13,9,0,5,6,12,11,10,7,8,7,2,8,5,12,2,5,8,13,6,7,4,5,12,10,12,7,4,6,5,7,10,8,6,14,7,9,5,8,6,9,7,13,11,10,10,11,15,10,3,8,9,1,
14,9,17,1,12,11,11,8,14,5,0,7,7,16,15,6,3,6,3,7,10,17,7,10,13,18,11,12,9,10,17,15,17,2,9,5,6,12,5,3,11,19,12,4,10,9,11,4,12,8,16,5,15,16,20,15,8,13,14,6,
15,10,10,6,7,4,8,15,15,6,7,0,14,17,16,1,6,9,4,2,11,18,4,7,6,19,6,5,8,11,18,16,18,7,4,8,7,13,12,10,12,20,13,3,3,10,12,5,11,15,17,12,16,9,21,10,9,10,9,5,
11,16,24,8,19,18,18,3,15,12,7,14,0,11,14,13,8,5,10,14,9,14,14,17,20,13,18,19,16,7,12,10,12,7,16,6,7,13,2,4,8,14,11,11,17,4,18,9,19,1,15,2,10,23,15,22,15,20,21,13,
4,9,17,17,12,15,11,8,8,11,16,17,11,0,7,18,19,10,13,19,6,7,13,10,13,4,11,12,9,6,1,1,5,18,13,11,10,6,11,13,5,3,4,20,14,7,11,20,12,10,8,13,3,16,4,15,8,13,14,12,
3,6,10,16,9,14,8,15,1,10,15,16,14,7,0,17,18,9,12,18,5,2,12,9,10,3,10,11,8,7,6,6,2,17,12,10,9,3,12,12,6,4,3,19,13,10,4,19,5,15,1,12,4,9,5,8,7,6,7,11,
16,11,11,5,8,5,9,14,16,7,6,1,13,18,17,0,5,8,5,1,12,19,5,8,7,20,7,6,9,12,19,17,19,6,5,7,8,14,11,9,13,21,14,2,4,11,13,4,12,14,18,11,17,10,22,11,10,11,10,6,
17,12,16,2,11,10,10,11,17,8,3,6,8,19,18,5,0,9,6,6,13,20,6,9,12,21,10,11,10,13,20,18,20,1,8,8,9,15,8,6,14,22,15,3,9,12,14,1,13,9,19,6,18,15,23,14,11,12,13,7,
8,11,19,7,14,13,13,6,10,7,6,9,5,10,9,8,9,0,5,9,4,11,9,12,15,12,13,14,11,4,11,9,11,8,11,1,2,8,3,3,5,13,6,10,12,3,13,10,14,6,10,3,9,18,14,17,10,15,16,8,
11,6,14,4,9,8,8,11,11,2,3,4,10,13,12,5,6,5,0,6,7,14,4,7,10,15,8,9,6,7,14,12,14,5,6,4,3,9,8,6,8,16,9,7,7,6,8,7,9,11,13,8,12,13,17,12,5,10,11,3,
17,12,10,6,9,4,10,15,17,8,7,2,14,19,18,1,6,9,6,0,13,20,6,9,8,21,8,7,10,13,20,18,20,7,6,8,9,15,12,10,14,22,15,3,5,12,14,5,13,15,19,12,18,11,23,12,11,12,11,7,
4,7,15,11,10,9,9,10,6,5,10,11,9,6,5,12,13,4,7,13,0,7,7,8,11,8,9,10,7,2,7,5,7,12,7,5,4,4,7,7,1,9,2,14,8,5,9,14,10,10,6,7,5,14,10,13,6,11,12,6,
3,8,10,18,11,16,10,15,3,12,17,18,14,7,2,19,20,11,14,20,7,0,14,11,12,3,12,13,10,7,6,6,2,19,14,12,11,5,12,14,6,4,5,21,15,10,6,21,7,15,1,14,4,9,3,8,9,8,9,13,
11,6,10,6,5,4,4,15,11,2,7,4,14,13,12,5,6,9,4,6,7,14,0,3,6,15,4,5,4,7,14,12,14,7,2,8,7,9,12,10,8,16,9,7,3,10,8,7,7,15,13,12,12,9,17,8,5,6,7,1,
8,3,7,9,2,5,1,18,8,5,10,7,17,10,9,8,9,12,7,9,8,11,3,0,3,12,1,2,1,10,11,9,11,10,3,11,10,6,15,13,9,13,6,10,4,13,5,10,4,18,10,15,9,6,14,5,2,3,4,4,
9,4,4,12,1,4,2,21,9,8,13,6,20,13,10,7,12,15,10,8,11,12,6,3,0,13,2,1,4,13,12,12,12,13,4,14,13,7,18,16,12,14,9,9,3,16,6,11,5,21,11,18,10,3,15,4,5,4,3,7,
4,9,13,19,12,17,11,12,4,13,18,19,13,4,3,20,21,12,15,21,8,3,15,12,13,0,13,14,11,8,3,3,1,20,15,13,12,6,13,15,7,1,6,22,16,9,7,22,8,12,4,15,3,12,2,11,10,9,10,14,
9,4,6,10,1,4,2,19,9,6,11,6,18,11,10,7,10,13,8,8,9,12,4,1,2,13,0,1,2,11,12,10,12,11,2,12,11,7,16,14,10,14,7,9,3,14,6,9,5,19,11,16,10,5,15,4,3,4,3,5,
10,5,5,11,2,3,3,20,10,7,12,5,19,12,11,6,11,14,9,7,10,13,5,2,1,14,1,0,3,12,13,11,13,12,3,13,12,8,17,15,11,15,8,8,2,15,7,10,6,20,12,17,11,4,16,5,4,5,4,6,
7,2,8,8,3,6,2,17,7,4,9,8,16,9,8,9,10,11,6,10,7,10,4,1,4,11,2,3,0,9,10,8,10,9,4,10,9,5,14,12,8,12,5,11,5,12,4,11,3,17,9,14,8,7,13,6,1,4,5,3,
4,9,17,11,12,11,11,8,8,5,10,11,7,6,7,12,13,4,7,13,2,7,7,10,13,8,11,12,9,0,7,5,7,12,9,5,4,6,5,7,1,9,4,14,10,3,11,14,12,8,8,7,5,16,10,15,8,13,14,6,
3,8,16,18,11,16,10,9,7,12,17,18,12,1,6,19,20,11,14,20,7,6,14,11,12,3,12,13,10,7,0,2,4,19,14,12,11,5,12,14,6,2,5,21,15,8,10,21,11,11,7,14,2,15,3,14,9,12,13,13,
3,8,16,16,11,14,10,9,7,10,15,16,10,1,6,17,18,9,12,18,5,6,12,9,12,3,10,11,8,5,2,0,4,17,12,10,9,5,10,12,4,4,3,19,13,6,10,19,11,9,7,12,2,15,5,14,7,12,13,11,
3,8,12,18,11,16,10,13,3,12,17,18,12,5,2,19,20,11,14,20,7,2,14,11,12,1,12,13,10,7,4,4,0,19,14,12,11,5,12,14,6,2,5,21,15,8,6,21,7,13,3,14,2,11,3,10,9,8,9,13,
16,11,17,1,12,11,11,10,16,7,2,7,7,18,17,6,1,8,5,7,12,19,7,10,13,20,11,12,9,12,19,17,19,0,9,7,8,14,7,5,13,21,14,4,10,11,13,2,12,8,18,5,17,16,22,15,10,13,14,6,
11,6,8,8,3,2,4,17,11,4,9,4,16,13,12,5,8,11,6,6,7,14,2,3,4,15,2,3,4,9,14,12,14,9,0,10,9,9,14,12,8,16,9,7,1,12,8,7,7,17,13,14,12,7,17,6,5,6,5,3,
9,10,18,6,13,12,12,7,9,6,5,8,6,11,10,7,8,1,4,8,5,12,8,11,14,13,12,13,10,5,12,10,12,7,10,0,1,7,4,2,6,14,7,9,11,4,12,9,13,7,11,4,10,17,15,16,9,14,15,7,
8,9,17,7,12,11,11,8,8,5,6,7,7,10,9,8,9,2,3,9,4,11,7,10,13,12,11,12,9,4,11,9,11,8,9,1,0,6,5,3,5,13,6,10,10,3,11,10,12,8,10,5,9,16,14,15,8,13,14,6,
2,3,11,13,6,11,5,14,2,7,12,13,13,6,3,14,15,8,9,15,4,5,9,6,7,6,7,8,5,6,5,5,5,14,9,7,6,0,11,9,5,7,2,16,10,9,5,16,6,14,4,11,3,10,8,9,4,7,8,8,
9,14,22,6,17,16,16,3,13,10,5,12,2,11,12,11,8,3,8,12,7,12,12,15,18,13,16,17,14,5,12,10,12,7,14,4,5,11,0,2,6,14,9,9,15,4,16,9,17,3,13,2,10,21,15,20,13,18,19,11,
11,12,20,4,15,14,14,5,11,8,3,10,4,13,12,9,6,3,6,10,7,14,10,13,16,15,14,15,12,7,14,12,14,5,12,2,3,9,2,0,8,16,9,7,13,6,14,7,15,5,13,2,12,19,17,18,11,16,17,9,
3,8,16,12,11,10,10,9,7,6,11,12,8,5,6,13,14,5,8,14,1,6,8,9,12,7,10,11,8,1,6,4,6,13,8,6,5,5,6,8,0,8,3,15,9,4,10,15,11,9,7,8,4,15,9,14,7,12,13,7,
5,10,14,20,13,18,12,11,5,14,19,20,14,3,4,21,22,13,16,22,9,4,16,13,14,1,14,15,12,9,2,4,2,21,16,14,13,7,14,16,8,0,7,23,17,10,8,23,9,13,5,16,4,13,1,12,11,10,11,15,
2,5,13,13,8,11,7,12,4,7,12,13,11,4,3,14,15,6,9,15,2,5,9,6,9,6,7,8,5,4,5,3,5,14,9,7,6,2,9,9,3,7,0,16,10,7,7,16,8,12,4,9,3,12,8,11,4,9,10,8,
18,13,13,3,10,7,11,12,18,9,4,3,11,20,19,2,3,10,7,3,14,21,7,10,9,22,9,8,11,14,21,19,21,4,7,9,10,16,9,7,15,23,16,0,6,13,15,2,14,12,20,9,19,12,24,13,12,13,12,8,
12,7,7,9,4,1,5,18,12,5,10,3,17,14,13,4,9,12,7,5,8,15,3,4,3,16,3,2,5,10,15,13,15,10,1,11,10,10,15,13,9,17,10,6,0,13,9,8,8,18,14,15,13,6,18,7,6,7,6,4,
7,12,20,10,15,14,14,5,11,8,9,10,4,7,10,11,12,3,6,12,5,10,10,13,16,9,14,15,12,3,8,6,8,11,12,4,3,9,4,6,4,10,7,13,13,0,14,13,15,5,11,6,6,19,11,18,11,16,17,9,
7,2,6,12,5,10,4,19,3,6,11,12,18,11,4,13,14,13,8,14,9,6,8,5,6,7,6,7,4,11,10,10,6,13,8,12,11,5,16,14,10,8,7,15,9,14,0,15,1,19,5,16,8,5,9,4,3,2,3,7,
18,13,15,3,10,9,11,12,18,9,4,5,9,20,19,4,1,10,7,5,14,21,7,10,11,22,9,10,11,14,21,19,21,2,7,9,10,16,9,7,15,23,16,2,8,13,15,0,14,10,20,7,19,14,24,13,12,13,12,8,
8,3,5,11,4,9,3,20,4,7,12,11,19,12,5,12,13,14,9,13,10,7,7,4,5,8,5,6,3,12,11,11,7,12,7,13,12,6,17,15,11,9,8,14,8,15,1,14,0,20,6,17,9,4,10,3,4,1,2,6,
12,17,25,9,20,19,19,2,16,13,8,15,1,10,15,14,9,6,11,15,10,15,15,18,21,12,19,20,17,8,11,9,13,8,17,7,8,14,3,5,9,13,12,12,18,5,19,10,20,0,16,3,11,24,14,23,16,21,22,14,
4,7,9,17,10,15,9,16,2,11,16,17,15,8,1,18,19,10,13,19,6,1,13,10,11,4,11,12,9,8,7,7,3,18,13,11,10,4,13,13,7,5,4,20,14,11,5,20,6,16,0,13,5,8,4,7,8,7,8,12,
11,14,22,6,17,16,16,5,13,10,5,12,2,13,12,11,6,3,8,12,7,14,12,15,18,15,16,17,14,7,14,12,14,5,14,4,5,11,2,2,8,16,9,9,15,6,16,7,17,3,13,0,12,21,17,20,13,18,19,11,
1,6,14,16,9,14,8,11,5,10,15,16,10,3,4,17,18,9,12,18,5,4,12,9,10,3,10,11,8,5,2,2,2,17,12,10,9,3,10,12,4,4,3,19,13,6,8,19,9,11,5,12,0,13,5,12,7,10,11,11,
12,7,1,15,4,7,5,24,8,11,16,9,23,16,9,10,15,18,13,11,14,9,9,6,3,12,5,4,7,16,15,15,11,16,7,17,16,10,21,19,15,13,12,12,6,19,5,14,4,24,8,21,13,0,12,1,8,3,2,10,
6,11,13,21,14,19,13,12,6,15,20,21,15,4,5,22,23,14,17,23,10,3,17,14,15,2,15,16,13,10,3,5,3,22,17,15,14,8,15,17,9,1,8,24,18,11,9,24,10,14,4,17,5,12,0,11,12,11,12,16,
11,6,2,14,3,8,4,23,7,10,15,10,22,15,8,11,14,17,12,12,13,8,8,5,4,11,4,5,6,15,14,14,10,15,6,16,15,9,20,18,14,12,11,13,7,18,4,13,3,23,7,20,12,1,11,0,7,2,1,9,
6,1,9,9,4,7,3,16,6,3,8,9,15,8,7,10,11,10,5,11,6,9,5,2,5,10,3,4,1,8,9,7,9,10,5,9,8,4,13,11,7,11,4,12,6,11,3,12,4,16,8,13,7,8,12,7,0,5,6,4,
9,4,4,12,3,8,2,21,5,8,13,10,20,13,6,11,12,15,10,12,11,8,6,3,4,9,4,5,4,13,12,12,8,13,6,14,13,7,18,16,12,10,9,13,7,16,2,13,1,21,7,18,10,3,11,2,5,0,1,7,
10,5,3,13,2,7,3,22,6,9,14,9,21,14,7,10,13,16,11,11,12,9,7,4,3,10,3,4,5,14,13,13,9,14,5,15,14,8,19,17,13,11,10,12,6,17,3,12,2,22,8,19,11,2,12,1,6,1,0,8,
10,5,11,5,6,5,5,14,10,1,6,5,13,12,11,6,7,8,3,7,6,13,1,4,7,14,5,6,3,6,13,11,13,6,3,7,6,8,11,9,7,15,8,8,4,9,7,8,6,14,12,11,11,10,16,9,4,7,8,0
	};
	static double D[ROBOT_NUM][TASK_NUM] =
	{
		16,11,3,19,8,9,9,28,12,15,20,13,27,20,13,14,19,22,17,13,18,13,13,10,7,16,9,8,11,20,19,19,15,20,11,21,20,14,25,23,19,17,16,16,10,23,9,18,8,28,12,25,17,4,16,5,12,7,6,14,
15,10,2,18,7,8,8,27,11,14,19,12,26,19,12,13,18,21,16,12,17,12,12,9,6,15,8,7,10,19,18,18,14,19,10,20,19,13,24,22,18,16,15,15,9,22,8,17,7,27,11,24,16,3,15,4,11,6,5,13,
14,9,3,17,6,9,7,26,10,13,18,11,25,18,11,12,17,20,15,13,16,11,11,8,5,14,7,6,9,18,17,17,13,18,9,19,18,12,23,21,17,15,14,14,8,21,7,16,6,26,10,23,15,2,14,3,10,5,4,12,
13,8,4,16,5,10,6,25,9,12,17,12,24,17,10,13,16,19,14,14,15,10,10,7,6,13,6,7,8,17,16,16,12,17,8,18,17,11,22,20,16,14,13,15,9,20,6,15,5,25,9,22,14,3,13,2,9,4,3,11,
12,7,5,15,6,11,5,24,8,11,16,13,23,16,9,14,15,18,13,15,14,9,9,6,7,12,7,8,7,16,15,15,11,16,9,17,16,10,21,19,15,13,12,16,10,19,5,16,4,24,8,21,13,4,12,3,8,3,4,10,
5,10,12,20,13,18,12,17,5,14,19,20,16,9,4,21,22,13,16,22,9,2,16,13,14,5,14,15,12,9,8,8,4,21,16,14,13,7,14,16,8,6,7,23,17,12,8,23,9,17,3,16,6,11,5,10,11,10,11,15,
6,11,13,21,14,19,13,16,6,15,20,21,15,8,5,22,23,14,17,23,10,3,17,14,15,4,15,16,13,10,7,7,3,22,17,15,14,8,15,17,9,5,8,24,18,11,9,24,10,16,4,17,5,12,4,11,12,11,12,16,
7,12,14,22,15,20,14,15,7,16,21,22,16,7,6,23,24,15,18,24,11,4,18,15,16,3,16,17,14,11,6,6,4,23,18,16,15,9,16,18,10,4,9,25,19,12,10,25,11,15,5,18,6,13,3,12,13,12,13,17,
8,13,15,23,16,21,15,14,8,17,22,23,17,6,7,24,25,16,19,25,12,5,19,16,17,4,17,18,15,12,5,7,5,24,19,17,16,10,17,19,11,3,10,26,20,13,11,26,12,16,6,19,7,14,2,13,14,13,14,18,
9,14,16,24,17,22,16,15,9,18,23,24,18,7,8,25,26,17,20,26,13,6,20,17,18,5,18,19,16,13,6,8,6,25,20,18,17,11,18,20,12,4,11,27,21,14,12,27,13,17,7,20,8,15,3,14,15,14,15,19,
17,12,20,4,15,14,14,11,17,8,3,10,8,19,18,9,4,9,6,10,13,20,10,13,16,21,14,15,12,13,20,18,20,3,12,8,9,15,8,6,14,22,15,7,13,12,14,5,15,9,19,6,18,19,23,18,11,16,17,9,
16,13,21,5,16,15,15,10,16,9,4,11,7,18,17,10,5,8,7,11,12,19,11,14,17,20,15,16,13,12,19,17,19,4,13,7,8,14,7,5,13,21,14,8,14,11,15,6,16,8,18,5,17,20,22,19,12,17,18,10,
15,14,22,6,17,16,16,9,15,10,5,12,6,17,16,11,6,7,8,12,11,18,12,15,18,19,16,17,14,11,18,16,18,5,14,6,7,13,6,4,12,20,13,9,15,10,16,7,17,7,17,4,16,21,21,20,13,18,19,11,
14,15,23,7,18,17,17,8,14,11,6,13,5,16,15,12,7,6,9,13,10,17,13,16,19,18,17,18,15,10,17,15,17,6,15,5,6,12,5,3,11,19,12,10,16,9,17,8,18,6,16,3,15,22,20,21,14,19,20,12,
13,16,24,8,19,18,18,7,15,12,7,14,4,15,14,13,8,5,10,14,9,16,14,17,20,17,18,19,16,9,16,14,16,7,16,6,7,13,4,4,10,18,11,11,17,8,18,9,19,5,15,2,14,23,19,22,15,20,21,13
	};
	static double h[TASK_NUM] = { 16,14,18,34,20,30,18,40,8,22,32,34,38,24,10,36,38,28,26,38,20,10,26,20,22,16,22,24,18,24,22,22,14,36,26,26,24,12,34,30,22,18,16,40,28,30,10,40,12,40,8,34,18,16,16,14,16,14,16,24 };
	double Cmax = BIG_M;  // Initialize the objective function value
	static double finishing_time[ANT_NUM][ROBOT_NUM];  // The completion time of each AGV
	static double final_finishing_time[ROBOT_NUM];  // The completion time of all AGVs corresponding to the optimal solution
	for (int i = 0; i != ROBOT_NUM; ++i)
	{
		final_finishing_time[i] = 0;
	}
	static double pheromones[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM];  // Pheromone matrix
	for (int i = 0; i != TASK_NUM * ROBOT_NUM; ++i)
	{
		for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
		{
			pheromones[i][j] = P0;  // Pheromone matrix initialization
		}
	}
	static ant_node final_solution[NODE_NUM];  // The global optimal solution found
	for (int i = 0;i != NODE_NUM; ++i)
	{
		final_solution[i].task = BIG_M;
	}
	vector<int>least_Cmax_combination_num;  // The binary number of the solution with the least Cmax is saved for global updates
	int current_ant_num = 0;
	int current_node_num = 0;
	int iteration_num = 0;
	int iteration_num_of_Cmax = 0;

	/*The core code of the ant colony algorithm*/
	for (int u = 0; u != T_MAX; ++u)
	{
		iteration_num = u;
		/*Before each iteration, reset the power of all the AGVs on the ants*/
		static int temp_M[ROBOT_NUM] = { 95,82,92,116,125,107,113,101,113,105,120,121,813,598,392 };  //BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M,BIG_M
		static double M[ANT_NUM][ROBOT_NUM];
		for (int a = 0; a != ANT_NUM; ++a)
		{
			for (int r = 0; r != ROBOT_NUM; ++r)
			{
				M[a][r] = temp_M[r];
			}
		}
		/*Before each iteration, reset the finish time of all ants*/
		for (int n = 0; n != ANT_NUM; ++n)
		{
			for (int m = 0; m != ROBOT_NUM; ++m)
			{
				finishing_time[n][m] = 0;
			}
		}

		/*Before each iteration, reset the node candidate set for all ants*/
		static ant_node candidate_set[ANT_NUM][TASK_NUM * ROBOT_NUM];
		for (int a = 0; a != ANT_NUM; ++a)
		{
			int combination_num = 0;
			for (int m = 0; m != TASK_NUM; ++m)
			{
				for (int n = 0; n != ROBOT_NUM; ++n)
				{
					candidate_set[a][combination_num].task = m;
					candidate_set[a][combination_num].AGV = n;
					candidate_set[a][combination_num].combination_num = combination_num;
					combination_num += 1;
				}
			}
		}

		vector<vector<vector<int>>>task_of_AGV;  // Record the task status of each AGV on the current ant, and clear it once for each iteration
		task_of_AGV.resize(ANT_NUM);
		for (int t = 0; t != ANT_NUM; ++t)
		{
			task_of_AGV[t].resize(ROBOT_NUM);
		}

		ant_node iteration_solution[ANT_NUM][NODE_NUM];  // Save the solutions found by all ants in one iteration
		vector<int>local_update_data_all_ant;  // Save the feasible binary numbers found by all ants in turn in each iteration
		vector<int>ant_num_of_local_update_data_all_ant;  // Save the ant numbers of all the viable binads found sequentially in each iteration
		vector<int>former_local_update_set;
		vector<int>ant_num_for_former_local_update_set;
		static int last_node_find_flag[ANT_NUM];
		for (int w = 0; w != ANT_NUM; ++w)
		{
			last_node_find_flag[w] = 1;  // Both initializations are 1, and when no viable binary is found, change to 0
		}
		for (int node_num=0;node_num!=NODE_NUM;++node_num)
		{
			current_node_num = node_num;
			vector<int>local_update_set;  // Saves the binary numbers of all ant selections for the current node for local updates
			vector<int>ant_num_for_local_update_set;  // Save the ant number where the current node finds a viable binad
			for (int ant_num = 0; ant_num != ANT_NUM; ++ant_num)  // All ants search the path in parallel
			{
				current_ant_num = ant_num;
				/*Assign a binary to the current node of the ant for which the previous node found a viable binad*/
				if (last_node_find_flag[current_ant_num] == 1)
				{
					ant_path_selection(iteration_solution, current_ant_num, current_node_num, candidate_set, pheromones, finishing_time, M, D, d, h, task_of_AGV, local_update_set, ant_num_for_local_update_set, last_node_find_flag);
				}
			}
			if (current_node_num != 0)  // Local updates start at the second node
			{
				/*The pheromones are updated locally only after the current node of all ants is determined*/
				local_pheromone_update(pheromones, local_update_set,former_local_update_set, last_node_find_flag,ant_num_for_former_local_update_set,ant_num_for_local_update_set);
				former_local_update_set.clear();
				ant_num_for_former_local_update_set.clear();
				for (int b = 0; b != local_update_set.size(); ++b)
				{
					former_local_update_set.push_back(local_update_set[b]);
					ant_num_for_former_local_update_set.push_back(ant_num_for_local_update_set[b]);
				}
			}
			else  // Save the binary numbers found by all ants on the first node
			{
				for (int b = 0; b != local_update_set.size(); ++b)
				{
					former_local_update_set.push_back(local_update_set[b]);
					ant_num_for_former_local_update_set.push_back(ant_num_for_local_update_set[b]);
				}
			}
			/*The binary index that needs to be updated for each ant at each step in the process is recorded for global update*/
			for (int m = 0; m != local_update_set.size(); ++m)
			{
				local_update_data_all_ant.push_back(local_update_set[m]);
				ant_num_of_local_update_data_all_ant.push_back(ant_num_for_local_update_set[m]);
			}
		}

		/*The global pheromones are updated at the end of each iteration*/
		global_pheromone_update(iteration_solution, final_solution, pheromones, local_update_data_all_ant, least_Cmax_combination_num, Cmax, finishing_time,task_of_AGV,iteration_num,iteration_num_of_Cmax,final_finishing_time,ant_num_of_local_update_data_all_ant);

	}
	finish = clock();
	double time_consuming = double(finish - start) / CLOCKS_PER_SEC;
	int finished_task_num = 0;
	for (int i = 0; i != NODE_NUM; ++i)
	{
		if (final_solution[i].task != BIG_M)
		{
			finished_task_num += 1;
		}
	}
	if (finished_task_num == TASK_NUM)
	{
		cout << "==============================================================" << endl;
		cout << "Successfully find a workable solution!" << endl;
		cout << "==============================================================" << endl;
		cout << "The solution time of the algorithm:" << time_consuming << "seconds" << endl;
		cout << "After the [ " << iteration_num_of_Cmax + 1 << " ] iteration, the global optimal Cmax = " << Cmax << " is found" << endl;
		cout << "Cmax corresponding to the global optimal solution: " << Cmax << endl;
		cout << "The completion time of each AGV in the global optimal solution is as follows: " << endl;
		for (int i = 0; i != ROBOT_NUM; ++i)
		{
			cout << "finishing_time[" << i + 1 << "] = " << final_finishing_time[i] << endl;
		}
		cout << "The final_solution of the global optimum: " << endl;
		for (int n = 0; n != NODE_NUM; ++n)
		{
			cout << "Node[" << n << "] : Binary " << final_solution[n].combination_num << ",（" << final_solution[n].task << "," << final_solution[n].AGV << "）" << endl;
		}
	}
	else
	{
		cout << "=============================================================" << endl;
		cout << "Hint: All iterations have been completed, but no workable solution has been found!" << endl;
		cout << "=============================================================" << endl;
	}
	getchar();
	return 0;
}

/*************************Function definition*************************/
//**Ant path selection function**//
void ant_path_selection(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM],int current_ant_num,int current_node_num,ant_node(&candidate_set)[ANT_NUM][TASK_NUM * ROBOT_NUM],double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM],double(&finishing_time)[ANT_NUM][ROBOT_NUM],double(&M)[ANT_NUM][ROBOT_NUM],double D[ROBOT_NUM][TASK_NUM],double d[TASK_NUM][TASK_NUM],double h[TASK_NUM], vector<vector<vector<int>>>(&task_of_AGV), vector<int>(&local_update_set), vector<int>(&ant_num_for_local_update_set),int(&last_node_find_flag)[ANT_NUM])
{
	/*Count the feasible node set of the current node of the current ant*/
	availiable_node availiable_set;

	if (current_node_num == 0)
	{
		for (int i = 0; i != TASK_NUM * ROBOT_NUM; ++i)
		{
			/*If it is the first node, all AGVs have not yet scheduled tasks*/
			int current_task = candidate_set[current_ant_num][i].task;
			int AGV_num = candidate_set[current_ant_num][i].AGV;
			if (M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task] > 0)
			{
				availiable_set.combination_num.push_back(candidate_set[current_ant_num][i].combination_num);
				availiable_set.task.push_back(candidate_set[current_ant_num][i].task);
				availiable_set.AGV.push_back(candidate_set[current_ant_num][i].AGV);
			}
		}
		if (availiable_set.combination_num.size() == 0)
		{
			/*If the current node of the current ant does not find a viable binad, the subsequent nodes do not need to be found*/
			last_node_find_flag[current_ant_num] = 0;
			return;
		}

		int random_num = random_int(0, availiable_set.combination_num.size() - 1);
		int next_combination_num = availiable_set.combination_num[random_num];
		int current_task = availiable_set.task[random_num];
		int AGV_num = availiable_set.AGV[random_num];
		iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
		iteration_solution[current_ant_num][current_node_num].task = current_task;
		iteration_solution[current_ant_num][current_node_num].AGV = AGV_num;
		finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];  // Update the completion time of the current AGV
		M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];  // Update the current power level of the AGV
		task_of_AGV[current_ant_num][AGV_num].push_back(current_task);  // Save the tasks assigned to the current AGV	
		local_update_set.push_back(next_combination_num);  // The index of the binary selected this time is recorded and used to update the pheromone locally
		ant_num_for_local_update_set.push_back(current_ant_num);
		last_node_find_flag[current_ant_num] = 1;
		for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
		{
			if (candidate_set[current_ant_num][j].task == current_task)
			{
				candidate_set[current_ant_num][j].combination_num = BIG_M;  // Tag the assigned binary
			}
		}
	}
	/*If it's not the ant's first node, you need to decide which method to use to choose the binary based on probability*/
	else
	{
		int temp_node_num = current_node_num - 1;
		int last_combination_num = iteration_solution[current_ant_num][temp_node_num].combination_num;
		int next_combination_num = 0;
		static double heuristic_information[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM];
		for (int p = 0; p != TASK_NUM * ROBOT_NUM; ++p)
		{
			for (int q = 0; q != TASK_NUM * ROBOT_NUM; ++q)
			{
				heuristic_information[p][q] = 0;
			}
		}
		double max_heuristic = 0;  // Maximum comprehensive heuristic information intensity
		double sum_heuristic = 0;
		int ant_Cmax_change = 0;
		double AGV_utilization_change = 0;
		for (int i = 0; i != TASK_NUM * ROBOT_NUM; ++i)
		{
			if (candidate_set[current_ant_num][i].combination_num != BIG_M)  // Exclude scheduled tasks
			{
				int current_task = candidate_set[current_ant_num][i].task;
				int AGV_num = candidate_set[current_ant_num][i].AGV;
				if (task_of_AGV[current_ant_num][AGV_num].size() == 0)  // If the current AGV has not scheduled a task
				{
					if (M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task] > 0)
					{
						availiable_set.combination_num.push_back(candidate_set[current_ant_num][i].combination_num);
						availiable_set.task.push_back(candidate_set[current_ant_num][i].task);
						availiable_set.AGV.push_back(candidate_set[current_ant_num][i].AGV);
					}
				}
				else  // If the current AGV has scheduled a task
				{
					int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
					int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
					if (M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task] > 0)
					{
						availiable_set.combination_num.push_back(candidate_set[current_ant_num][i].combination_num);
						availiable_set.task.push_back(candidate_set[current_ant_num][i].task);
						availiable_set.AGV.push_back(candidate_set[current_ant_num][i].AGV);
					}
				}
			}
		}
		if (availiable_set.combination_num.size() == 0)
		{
			/*If the current node of the current ant does not find a viable binary, the subsequent nodes do not need to be found*/
			last_node_find_flag[current_ant_num] = 0;
			return;
		}
		
		for (int i = 0; i != availiable_set.combination_num.size(); ++i)
		{
			int temp_next_combination_num = availiable_set.combination_num[i];
			int current_task = availiable_set.task[i];
			int current_AGV = availiable_set.AGV[i];
			/*Calculate the Cmax change of the current AGV caused by the current binary*/
			ant_Cmax_change = Cmax_change_statistics(iteration_solution, finishing_time, current_ant_num, current_node_num, current_AGV, current_task, D, d, h);
			/*Calculate the change in utilization of the current AGV caused by the current binary*/
			AGV_utilization_change = AGV_utilization_change_statistics(iteration_solution, current_ant_num, current_node_num, current_AGV, current_task, D, d, h);
			heuristic_information[last_combination_num][temp_next_combination_num] = pow(pheromones[last_combination_num][temp_next_combination_num], A) * pow(1 / (ant_Cmax_change + 0.5), B) * exp(AGV_utilization_change * C);
			sum_heuristic += heuristic_information[last_combination_num][temp_next_combination_num]; 
			if (heuristic_information[last_combination_num][temp_next_combination_num] > max_heuristic)
			{
				max_heuristic = heuristic_information[last_combination_num][temp_next_combination_num];
			}
		}

		double check_num = random_double(0, 1);
		if (check_num > Q0)  // If the random probability is greater than q0, the next binary is randomly selected according to the method of roulette
		{
			vector<int>dij_combination_num;
			vector<double>dij_rate;
			vector<double>dij_cumulate_rate;
			double temp_cumulate_rate = 0;
			for (int i = 0; i != availiable_set.combination_num.size(); ++i)
			{
				int temp_next_combination_num = availiable_set.combination_num[i];
				dij_combination_num.push_back(temp_next_combination_num);
				double temp_rate = heuristic_information[last_combination_num][temp_next_combination_num] / sum_heuristic;
				dij_rate.push_back(temp_rate);
				temp_cumulate_rate = temp_cumulate_rate + temp_rate;
				dij_cumulate_rate.push_back(temp_cumulate_rate);
			}

			/*If there is only one viable binary*/
			if (availiable_set.combination_num.size() == 1)
			{
				/*Important: The binary number you want to select for the next node has been found*/
				next_combination_num = dij_combination_num[0];
				iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
				iteration_solution[current_ant_num][current_node_num].task = availiable_set.task[0];
				iteration_solution[current_ant_num][current_node_num].AGV = availiable_set.AGV[0];
				int current_task = iteration_solution[current_ant_num][current_node_num].task;
				int AGV_num = iteration_solution[current_ant_num][current_node_num].AGV;
				/*Check whether the AGV of the previous node is the same as that of the current node*/
				int last_node_num = current_node_num - 1;
				if (AGV_num == iteration_solution[current_ant_num][last_node_num].AGV)
				{
					/*Other tasks have been scheduled for the current AGV*/
					int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
					int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
					finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];  // Update the completion time of the current AGV
					M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];  // Update the current power level of the AGV
					task_of_AGV[current_ant_num][AGV_num].push_back(current_task);  // Save the tasks assigned to the current AGV
					local_update_set.push_back(next_combination_num);  // The number of the binary selected this time is recorded and used to update the pheromone locally
					ant_num_for_local_update_set.push_back(current_ant_num);
					last_node_find_flag[current_ant_num] = 1;
					for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
					{
						if (candidate_set[current_ant_num][j].task == current_task)
						{
							candidate_set[current_ant_num][j].combination_num = BIG_M;  // Tag the assigned binary
						}
					}
				}
				else  // The AGV of the current node is different from the AGV of the previous node
				{
					if (task_of_AGV[current_ant_num][AGV_num].size() == 0)
					{
						finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];
						M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];  // Update the current power level of the AGV
						task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
						local_update_set.push_back(next_combination_num);
						ant_num_for_local_update_set.push_back(current_ant_num);
						last_node_find_flag[current_ant_num] = 1;
						for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
						{
							if (candidate_set[current_ant_num][j].task == current_task)
							{
								candidate_set[current_ant_num][j].combination_num = BIG_M;
							}
						}
					}
					else
					{
						int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
						int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
						finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
						M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
						task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
						local_update_set.push_back(next_combination_num);
						ant_num_for_local_update_set.push_back(current_ant_num);
						last_node_find_flag[current_ant_num] = 1;
						for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
						{
							if (candidate_set[current_ant_num][j].task == current_task)
							{
								candidate_set[current_ant_num][j].combination_num = BIG_M;
							}
						}
					}
				}
			}
			/*If there are more than one viable binary*/
			else
			{
				/*By roulette, the choice is made randomly based on the probability of the strength of the combined heuristic information of each binary group*/
				double temp_random_num = random_double(0, 1);
				for (int i = 0; i != dij_cumulate_rate.size() - 1; ++i)
				{
					if (temp_random_num <= dij_cumulate_rate[0])
					{
						next_combination_num = dij_combination_num[0];
						iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
						iteration_solution[current_ant_num][current_node_num].task = availiable_set.task[0];
						iteration_solution[current_ant_num][current_node_num].AGV = availiable_set.AGV[0];
						int current_task = iteration_solution[current_ant_num][current_node_num].task;
						int AGV_num = iteration_solution[current_ant_num][current_node_num].AGV;
						/*Check whether the AGV of the previous node is the same as that of the current node*/
						int last_node_num = current_node_num - 1;
						if (AGV_num == iteration_solution[current_ant_num][last_node_num].AGV)
						{
							/*Other tasks have been scheduled for the current AGV*/
							int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
							int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
							finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
							M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
							task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
							local_update_set.push_back(next_combination_num);
							ant_num_for_local_update_set.push_back(current_ant_num);
							last_node_find_flag[current_ant_num] = 1;
							for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
							{
								if (candidate_set[current_ant_num][j].task == current_task)
								{
									candidate_set[current_ant_num][j].combination_num = BIG_M;
								}
							}
							break;  // As soon as the binary of the current node of the current ant is determined, the for loop jumps out
						}
						else
						{
							if (task_of_AGV[current_ant_num][AGV_num].size() == 0) 
							{
								finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];
								M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];
								task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
								local_update_set.push_back(next_combination_num);
								ant_num_for_local_update_set.push_back(current_ant_num);
								last_node_find_flag[current_ant_num] = 1;
								for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
								{
									if (candidate_set[current_ant_num][j].task == current_task)
									{
										candidate_set[current_ant_num][j].combination_num = BIG_M;
									}
								}
								break;
							}
							else
							{
								int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
								int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
								finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
								M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
								task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
								local_update_set.push_back(next_combination_num);
								ant_num_for_local_update_set.push_back(current_ant_num);
								last_node_find_flag[current_ant_num] = 1;
								for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
								{
									if (candidate_set[current_ant_num][j].task == current_task)
									{
										candidate_set[current_ant_num][j].combination_num = BIG_M;
									}
								}
								break;
							}
						}
					}
					else
					{
						int temp_num_dij = i + 1;
						if (( temp_random_num > dij_cumulate_rate[i]) && (temp_random_num <= dij_cumulate_rate[temp_num_dij]))
						{
							next_combination_num = dij_combination_num[temp_num_dij];  // Select the i+1 binary
							iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
							iteration_solution[current_ant_num][current_node_num].task = availiable_set.task[temp_num_dij];
							iteration_solution[current_ant_num][current_node_num].AGV = availiable_set.AGV[temp_num_dij];
							int current_task = iteration_solution[current_ant_num][current_node_num].task;
							int AGV_num = iteration_solution[current_ant_num][current_node_num].AGV;
							int last_node_num = current_node_num - 1;
							if (AGV_num == iteration_solution[current_ant_num][last_node_num].AGV)
							{
								int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
								int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
								finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
								M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
								task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
								local_update_set.push_back(next_combination_num);
								ant_num_for_local_update_set.push_back(current_ant_num);
								last_node_find_flag[current_ant_num] = 1;
								for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
								{
									if (candidate_set[current_ant_num][j].task == current_task)
									{
										candidate_set[current_ant_num][j].combination_num = BIG_M;
									}
								}
								break;
							}
							else
							{
								if (task_of_AGV[current_ant_num][AGV_num].size() == 0)
								{
									finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];
									M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];
									task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
									local_update_set.push_back(next_combination_num);
									ant_num_for_local_update_set.push_back(current_ant_num);
									last_node_find_flag[current_ant_num] = 1;
									for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
									{
										if (candidate_set[current_ant_num][j].task == current_task)
										{
											candidate_set[current_ant_num][j].combination_num = BIG_M;
										}
									}
									break;
								}
								else
								{
									int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
									int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
									finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
									M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
									task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
									local_update_set.push_back(next_combination_num);
									ant_num_for_local_update_set.push_back(current_ant_num);
									last_node_find_flag[current_ant_num] = 1;
									for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
									{
										if (candidate_set[current_ant_num][j].task == current_task)
										{
											candidate_set[current_ant_num][j].combination_num = BIG_M;
										}
									}
									break;
								}
							}
						}
						else
						{
							if (temp_random_num >= 1)
							{
								int temp_num_last_position = dij_combination_num.size() - 1;
								next_combination_num = dij_combination_num[temp_num_last_position];  //选择最后一个二元组
								iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
								iteration_solution[current_ant_num][current_node_num].task = availiable_set.task[temp_num_last_position];
								iteration_solution[current_ant_num][current_node_num].AGV = availiable_set.AGV[temp_num_last_position];
								int current_task = iteration_solution[current_ant_num][current_node_num].task;
								int AGV_num = iteration_solution[current_ant_num][current_node_num].AGV;
								int last_node_num = current_node_num - 1;
								if (AGV_num == iteration_solution[current_ant_num][last_node_num].AGV)
								{
									int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
									int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
									finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
									M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
									task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
									local_update_set.push_back(next_combination_num);
									ant_num_for_local_update_set.push_back(current_ant_num);
									last_node_find_flag[current_ant_num] = 1;
									for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
									{
										if (candidate_set[current_ant_num][j].task == current_task)
										{
											candidate_set[current_ant_num][j].combination_num = BIG_M;
										}
									}
									break;
								}
								else
								{
									if (task_of_AGV[current_ant_num][AGV_num].size() == 0)
									{
										finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];
										M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];
										task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
										local_update_set.push_back(next_combination_num);
										ant_num_for_local_update_set.push_back(current_ant_num);
										last_node_find_flag[current_ant_num] = 1;
										for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
										{
											if (candidate_set[current_ant_num][j].task == current_task)
											{
												candidate_set[current_ant_num][j].combination_num = BIG_M;
											}
										}
										break;
									}
									else
									{
										int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
										int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
										finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
										M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
										task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
										local_update_set.push_back(next_combination_num);
										ant_num_for_local_update_set.push_back(current_ant_num);
										last_node_find_flag[current_ant_num] = 1;
										for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
										{
											if (candidate_set[current_ant_num][j].task == current_task)
											{
												candidate_set[current_ant_num][j].combination_num = BIG_M;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
		else
		{
			for (int i = 0; i != availiable_set.combination_num.size(); ++i)
			{
				int temp_next_combination_num = availiable_set.combination_num[i];
				if (heuristic_information[last_combination_num][temp_next_combination_num] == max_heuristic)
				{
					next_combination_num = temp_next_combination_num;
					iteration_solution[current_ant_num][current_node_num].combination_num = next_combination_num;
					iteration_solution[current_ant_num][current_node_num].task = availiable_set.task[i];
					iteration_solution[current_ant_num][current_node_num].AGV = availiable_set.AGV[i];
					int current_task = iteration_solution[current_ant_num][current_node_num].task;
					int AGV_num = iteration_solution[current_ant_num][current_node_num].AGV;
					int last_node_num = current_node_num - 1;
					if (AGV_num == iteration_solution[current_ant_num][last_node_num].AGV)
					{
						int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
						int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
						finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
						M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
						task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
						local_update_set.push_back(next_combination_num);
						ant_num_for_local_update_set.push_back(current_ant_num);
						last_node_find_flag[current_ant_num] = 1;
						for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
						{
							if (candidate_set[current_ant_num][j].task == current_task)
							{
								candidate_set[current_ant_num][j].combination_num = BIG_M;
							}
						}
						break;
					}
					else
					{
						if (task_of_AGV[current_ant_num][AGV_num].size() == 0)
						{
							finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + D[AGV_num][current_task] + h[current_task];
							M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - D[AGV_num][current_task] - h[current_task];
							task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
							local_update_set.push_back(next_combination_num);
							ant_num_for_local_update_set.push_back(current_ant_num);
							last_node_find_flag[current_ant_num] = 1;
							for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
							{
								if (candidate_set[current_ant_num][j].task == current_task)
								{
									candidate_set[current_ant_num][j].combination_num = BIG_M;
								}
							}
							break;
						}
						else
						{
							int temp_num = task_of_AGV[current_ant_num][AGV_num].size() - 1;
							int last_task = task_of_AGV[current_ant_num][AGV_num][temp_num];
							finishing_time[current_ant_num][AGV_num] = finishing_time[current_ant_num][AGV_num] + d[last_task][current_task] + h[current_task];
							M[current_ant_num][AGV_num] = M[current_ant_num][AGV_num] - d[last_task][current_task] - h[current_task];
							task_of_AGV[current_ant_num][AGV_num].push_back(current_task);
							local_update_set.push_back(next_combination_num);
							ant_num_for_local_update_set.push_back(current_ant_num);
							last_node_find_flag[current_ant_num] = 1;
							for (int j = 0; j != TASK_NUM * ROBOT_NUM; ++j)
							{
								if (candidate_set[current_ant_num][j].task == current_task)
								{
									candidate_set[current_ant_num][j].combination_num = BIG_M;
								}
							}
							break;
						}
					}
				}
			}
		}
	}
}

//**Local pheromone update function**//
void local_pheromone_update(double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM], vector<int>(&local_update_set), vector<int>former_local_update_set, int(&last_node_find_flag)[ANT_NUM],vector<int>(&ant_num_for_former_local_update_set), vector<int>(&ant_num_for_local_update_set))
{
	int finishing_ant_num = 0;
	/*It is only necessary to count the number of ants that the previous node found a viable binary and the number of ants that the current node also found a binad*/
	for (int i = 0; i != ant_num_for_former_local_update_set.size(); ++i)
	{
		if (last_node_find_flag[ant_num_for_former_local_update_set[i]] == 1)
		{
			finishing_ant_num = finishing_ant_num + 1;
		}
	}
	if (finishing_ant_num == ant_num_for_former_local_update_set.size())
	{
		for (int i = 0; i != local_update_set.size(); ++i)
		{
			int last_combination_num = former_local_update_set[i];
			int next_combination_num = local_update_set[i];
			pheromones[last_combination_num][next_combination_num] = (1 - R_L) * pheromones[last_combination_num][next_combination_num] + R_L * P0;  // Local pheromone updates

		}
	}
	else
	{
		vector<int>new_former_local_update_set;
		for (int i = 0; i != former_local_update_set.size(); ++i)
		{
			if (last_node_find_flag[ant_num_for_former_local_update_set[i]] == 1)
			{
				/*Saves the binary numbers found by the previous node for all ants found this time for local updates*/
				new_former_local_update_set.push_back(former_local_update_set[i]);
			}
		}
		for (int i = 0; i != local_update_set.size(); ++i)
		{
			int last_combination_num = new_former_local_update_set[i];
			int next_combination_num = local_update_set[i];
			pheromones[last_combination_num][next_combination_num] = (1 - R_L) * pheromones[last_combination_num][next_combination_num] + R_L * P0;  // Local pheromone updates
		}
	}
}

//**Global pheromone update function**//
void global_pheromone_update(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM],ant_node(&final_solution)[NODE_NUM],double(&pheromones)[TASK_NUM * ROBOT_NUM][TASK_NUM * ROBOT_NUM],vector<int>(&local_update_data_all_ant),vector<int>(&least_Cmax_combination_num),double(&Cmax),double(&finishing_time)[ANT_NUM][ROBOT_NUM], vector<vector<vector<int>>>task_of_AGV, int iteration_num, int(&iteration_num_of_Cmax), double(&final_finishing_time)[ROBOT_NUM], vector<int>(&ant_num_of_local_update_data_all_ant))
{
	vector<int>local_Cmax_set;
	vector<int>local_Cmax_ant_num;
	static int task_sum[ANT_NUM];
	for (int t = 0; t != ANT_NUM; ++t)
	{
		task_sum[t] = 0;
	}
	for (int a = 0; a != ANT_NUM; ++a)
	{
		for (int u = 0; u != ROBOT_NUM; ++u)
		{
			task_sum[a] = task_sum[a] + task_of_AGV[a][u].size();
		}
	}

	for (int a = 0; a != ANT_NUM; ++a)
	{
		/*Only the Cmax of all ants that have completed all node searches are counted*/
		if (task_sum[a] == TASK_NUM)  // Exclude the effect of Cmax in ants that did not complete the search
		{
			int current_ant_Cmax = 0;
			for (int r = 0; r != ROBOT_NUM; ++r)
			{
				if (finishing_time[a][r] > current_ant_Cmax) 
				{
					current_ant_Cmax = finishing_time[a][r];
				}
			}
			local_Cmax_set.push_back(current_ant_Cmax);
			local_Cmax_ant_num.push_back(a);
		}
	}
	if (local_Cmax_set.size() == 0)  // None of the ants found a complete workable solution
	{
		if (least_Cmax_combination_num.size() != 0)
		{
			for (int i = 0; i != least_Cmax_combination_num.size() - 1; ++i)
			{
				int last_combination_num = least_Cmax_combination_num[i];
				int temp_num = i + 1;
				int next_combination_num = least_Cmax_combination_num[temp_num];
				/*Global pheromone updates*/
				pheromones[last_combination_num][next_combination_num] = (1 - R_G) * pheromones[last_combination_num][next_combination_num] + R_G * (1 / Cmax);
			}
		}
	}
	else  // Ants find a complete workable solution
	{
		int temp_global_Cmax = BIG_M;
		int temp_global_Cmax_ant_num = 0;
		for (int i = 0; i != local_Cmax_set.size(); ++i)
		{
			if (local_Cmax_set[i] < temp_global_Cmax)
			{
				temp_global_Cmax = local_Cmax_set[i];
				temp_global_Cmax_ant_num = local_Cmax_ant_num[i];
			}
		}

		/*If the minimum Cmax found in this iteration is smaller than the previously found, it is saved as the global minimum Cmax and the pheromones are updated*/
		if (temp_global_Cmax < Cmax)
		{
			iteration_num_of_Cmax = iteration_num;
			Cmax = temp_global_Cmax;
			for (int i = 0; i != ROBOT_NUM; ++i)
			{
				final_finishing_time[i] = finishing_time[temp_global_Cmax_ant_num][i];
			}
			vector<int>global_update_set;
			for (int i = 0; i != local_update_data_all_ant.size(); ++i)
			{
				if (ant_num_of_local_update_data_all_ant[i] == temp_global_Cmax_ant_num)
				{
					global_update_set.push_back(local_update_data_all_ant[i]);
				}
			}
			least_Cmax_combination_num.clear();
			for (int i = 0; i != global_update_set.size(); ++i)
			{
				/*Update the binary index corresponding to the global optimal Cmax*/
				least_Cmax_combination_num.push_back(global_update_set[i]);
			}
			for (int i = 0; i != least_Cmax_combination_num.size() - 1; ++i)
			{
				int last_combination_num = least_Cmax_combination_num[i];
				int temp_num = i + 1;
				int next_combination_num = least_Cmax_combination_num[temp_num];
				/*Global pheromone updates*/
				pheromones[last_combination_num][next_combination_num] = (1 - R_G) * pheromones[last_combination_num][next_combination_num] + R_G * (1 / Cmax);
			}
			for (int i = 0; i != NODE_NUM; ++i)
			{
				final_solution[i].combination_num = iteration_solution[temp_global_Cmax_ant_num][i].combination_num;
				final_solution[i].task = iteration_solution[temp_global_Cmax_ant_num][i].task;
				final_solution[i].AGV = iteration_solution[temp_global_Cmax_ant_num][i].AGV;
			}
		}
		/*If the minimum Cmax found in this iteration is not as small as the previously saved global minimum Cmax, the pheromone corresponding to the original global minimum Cmax is still updated*/
		else
		{
			for (int i = 0; i != least_Cmax_combination_num.size() - 1; ++i)
			{
				int last_combination_num = least_Cmax_combination_num[i];
				int temp_num = i + 1;
				int next_combination_num = least_Cmax_combination_num[temp_num];
				/*Global pheromone updates*/
				pheromones[last_combination_num][next_combination_num] = (1 - R_G) * pheromones[last_combination_num][next_combination_num] + R_G * (1 / Cmax);
			}
		}
	}
	
}

//**Cmax change statistical function**//
double Cmax_change_statistics(ant_node(&solution)[ANT_NUM][NODE_NUM],double(&finishing_time)[ANT_NUM][ROBOT_NUM],int current_ant_num, int current_node_num, int current_AGV, int current_task, double D[ROBOT_NUM][TASK_NUM], double d[TASK_NUM][TASK_NUM], double h[TASK_NUM])
{
	static int temp_finishing_time[ROBOT_NUM];
	for (int r = 0; r != ROBOT_NUM; ++r)
	{
		temp_finishing_time[r] = finishing_time[current_ant_num][r];
	}
	for (int r = 0; r != ROBOT_NUM; ++r)
	{
		if (r==current_AGV)  // Updates the completion time of the AGV corresponding to the current binary
		{
			vector<int>task_of_current_AGV;  // Save the scheduled tasks on the current AGV
			for (int i = 0; i != current_node_num; ++i)
			{
				if (solution[current_ant_num][i].AGV == current_AGV)  // Find all the tasks that the current ant has scheduled on the current AGV
				{
					task_of_current_AGV.push_back(solution[current_ant_num][i].task);
				}
			}
			if (task_of_current_AGV.size() == 0)
			{
				temp_finishing_time[current_AGV] = D[current_AGV][current_task] + h[current_task] + temp_finishing_time[current_AGV];
			}
			else
			{
				int temp_num = task_of_current_AGV.size() - 1;
				temp_finishing_time[current_AGV] = d[task_of_current_AGV[temp_num]][current_task] + h[current_task] + temp_finishing_time[current_AGV];
			}
		}
	}
	int old_Cmax = 0;
	int new_Cmax = 0;
	int Cmax_change = 0;
	for (int r = 0; r != ROBOT_NUM; ++r)
	{
		if (finishing_time[current_ant_num][r] > old_Cmax)
		{
			old_Cmax = finishing_time[current_ant_num][r];
		}
	}
	for (int r = 0; r != ROBOT_NUM; ++r)
	{
		if (temp_finishing_time[r] > new_Cmax)
		{
			new_Cmax = temp_finishing_time[r];
		}
	}
	Cmax_change = new_Cmax - old_Cmax;
	return Cmax_change;
}

//**Statistical function of AGV utilization rate change**//
double AGV_utilization_change_statistics(ant_node(&iteration_solution)[ANT_NUM][NODE_NUM], int current_ant_num, int current_node_num, int current_AGV, int current_task, double D[ROBOT_NUM][TASK_NUM],double d[TASK_NUM][TASK_NUM], double h[TASK_NUM])
{
	double old_utilization = 0;
	double new_utilization = 0;
	double utilization_change = 0;
	vector<int>task_of_current_AGV;  // Save the scheduled tasks on the current AGV
	for (int i = 0; i != current_node_num; ++i)
	{
		if (iteration_solution[current_ant_num][i].AGV == current_AGV)  // Find all the tasks that the current ant has scheduled on the current AGV
		{
			task_of_current_AGV.push_back(iteration_solution[current_ant_num][i].task);
		}
	}
	
	if (task_of_current_AGV.size() == 0)
	{
		double molecule_1 = h[current_task]; 
		double denominator_1 = D[current_AGV][current_task] + h[current_task]; 
		new_utilization = molecule_1 / denominator_1;
		utilization_change = new_utilization;
	}
	else
	{
		double sum_handling_distance = 0;
		for (int k = 0; k != task_of_current_AGV.size(); ++k)  // Calculate the handling rack mileage of the current AGV
		{
			sum_handling_distance = sum_handling_distance + h[task_of_current_AGV[k]];
		}
		double distance = 0;
		/*If the current AGV has only scheduled one task before, distance=0*/
		for (int k = 1; k != task_of_current_AGV.size(); ++k)  // Calculate the total mileage of the current AGV except for the first shelf
		{
			int temp_num = k - 1;
			distance = distance + d[task_of_current_AGV[temp_num]][task_of_current_AGV[k]] + h[task_of_current_AGV[k]];
		}

		double old_molecule = sum_handling_distance;
		double old_denominator = D[current_AGV][task_of_current_AGV[0]] + h[task_of_current_AGV[0]] + distance;
		old_utilization = old_molecule / old_denominator;
		double new_molecule = sum_handling_distance + h[current_task];
		int temp_num = task_of_current_AGV.size() - 1;
		double new_denominator = D[current_AGV][task_of_current_AGV[0]] + h[task_of_current_AGV[0]] + distance + d[task_of_current_AGV[temp_num]][current_task] + h[current_task]; 
		new_utilization = new_molecule / new_denominator;
		utilization_change = new_utilization - old_utilization;
	}
	return utilization_change;
}