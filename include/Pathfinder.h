#pragma once
#ifndef PATHFINDER_H
#define PATHFINDER_H
#include "math.h"
#include <vector>
#include <chrono>
#include <cstring>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <tuple>
#include <utility>
#include<array>
using namespace std;

class Pathfinder
{
public:
	Pathfinder()
	{

	};
	~Pathfinder()
	{

	};
	// Creating a shortcut for int, int pair type
	typedef pair<int, int> Pair;
	// Creating a shortcut for tuple<int, int, int> type
	typedef tuple<double, int, int> Tuple;


	
	void help(const vector<vector<int>>& grid, const Pair& src, const Pair& dest)
	{
		std::cout << "hallo " << grid[0][0] << " " << src.first << " " << dest.second << std::endl;
	}
	// A Function to find the shortest path between a given
	// source cell to a destination cell according to A* Search
	// Algorithm
	
	void aStarSearch(const vector<vector<int>>& grid, const Pair& src, const Pair& dest, vector<std::pair<int, int>>& res)
	{
		//std::cout << "aStarSearch" << std::endl;
		ROW = grid.size();
		COL = grid[0].size();

		std::array<Pair, 3> eastPosib = { Pair(1,0),Pair(1,1),Pair(1,-1) };
		std::array<Pair, 3> eastPosib_start = { Pair(1,0),Pair(0,1),Pair(0,-1) };

		//zurück
		std::array<Pair, 3> westPosib = { Pair(-1,0),Pair(-1,1),Pair(-1,-1) };
		std::array<Pair, 3> westPosib_start = { Pair(1,0),Pair(0,1),Pair(0,-1) };


		std::array<Pair, 3> southPosib = { Pair(0,1),Pair(-1,1),Pair(+1,1) };
		std::array<Pair, 3> southPosib_start = { Pair(0,1),Pair(-1,0),Pair(+1,0) };

		//zurück
		std::array<Pair, 3> northPosib = { Pair(0,-1),Pair(-1,-1),Pair(+1,-1) };
		std::array<Pair, 3> northPosib_start = { Pair(0,1),Pair(-1,0),Pair(+1,0) };

		direction dir;
		if (abs(src.first - dest.first) > abs(src.second - dest.second))//wo ist die hauptrichtung
		{
			if (src.first < dest.first)
			{
				dir = direction::east;
				//std::cout << "links nach rechts east" << std::endl;
			}
			else
			{
				dir = direction::west;
				//std::cout << "rechts nach links west" << std::endl;
			}
				
		}
		else
		{
			if (src.second < dest.second)
			{
				dir = direction::south;
				//std::cout << "oben nach unten sout" << std::endl;
			}
			else
			{
				dir = direction::north;
			//	std::cout << "unten nach oben north" << std::endl;
			}
				
		}
		//std::cout << "direction: " << dir << std::endl;
		
		
		
		// If the source is out of range
		if (!isValid(grid, src)) {
			printf("Source is invalid\n");
			return;
		}

		// If the destination is out of range
		if (!isValid(grid, dest)) {
			printf("Destination is invalid\n");
			return;
		}

		// If the destination cell is the same as source cell
		if (isDestination(src, dest)) {
			printf("We are already at the destination\n");
			return;
		}
		// Either the source or the destination is blocked
		if (!isUnBlocked(grid, src)
			|| !isUnBlocked(grid, dest)) {
			
			printf("Source or the destination is blocked\n");
			return;
		}
		// Create a closed list and initialise it to false which
		// means that no cell has been included yet This closed
		// list is implemented as a boolean 2D vector
		std::vector<std::vector<bool>> closedList;// [ROW] [COL] ;
		closedList.resize(ROW);
		for (int r = 0; r < ROW; r++)
		{
			closedList[r].resize(COL);
			memset(&closedList[r][0], false, closedList[r].size() * sizeof closedList[r][0]);
		}
			



		//memset(closedList, false, sizeof(closedList));

		// Declare a 2D vector of structure to hold the details
		// of that cell
		vector<vector<cell>> cellDetails;

		cellDetails.resize(ROW);
		for (int r = 0; r < ROW; r++)
			cellDetails[r].resize(COL);

		int i, j;
		// Initialising the parameters of the starting node
		i = src.first, j = src.second;
		cellDetails[i][j].f = 0.0;
		cellDetails[i][j].g = 0.0;
		cellDetails[i][j].h = 0.0;
		cellDetails[i][j].parent = { i, j };

		/*
		Create an open list having information as-
		<f, <i, j>>
		where f = g + h,
		and i, j are the row and column index of that cell
		Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
		This open list is implemented as a set of tuple.*/
		std::priority_queue<Tuple, std::vector<Tuple>,
			std::greater<Tuple> >
			openList;

		// Put the starting cell on the open list and set its
		// 'f' as 0
		openList.emplace(0.0, i, j);

		// We set this boolean value as false as initially
		// the destination is not reached.
		while (!openList.empty()) {
			const Tuple& p = openList.top();
			// Add this vertex to the closed list
			i = get<1>(p); // second element of tupla
			j = get<2>(p); // third element of tupla
			Pair position = { i, j };
			// Remove this vertex from the open list
			openList.pop();
			closedList[i][j] = true;
			
			/*if (position == back)
			{
				if (dir == direction::north)
					dir = direction::south;

				if (dir == direction::north)
					dir = direction::south;

				if (dir == direction::west)
					dir = direction::east;

				if (dir == direction::east)
					dir = direction::west;
			}*/
			/*
					Generating all the 8 successor of this cell
							N.W N N.E
							\ | /
							\ | /
							W----Cell----E
									/ | \
							/ | \
							S.W S S.E

					Cell-->Popped Cell (i, j)
					N --> North	 (i-1, j)
					S --> South	 (i+1, j)
					E --> East	 (i, j+1)
					W --> West		 (i, j-1)
					N.E--> North-East (i-1, j+1)
					N.W--> North-West (i-1, j-1)
					S.E--> South-East (i+1, j+1)
					S.W--> South-West (i+1, j-1)
			*/
			//for (int add_x = -1; add_x <= 1; add_x++) {
			//	for (int add_y = -1; add_y <= 1; add_y++) {
			for(int l = 0; l < 3; l++ )
			{
				
				Pair neighbour;
				if (dir == direction::south)
				{
					neighbour = {i+ southPosib[l].first,j+southPosib[l].second };
					//start variable
					if (j == src.second)
					{
						neighbour = { i + southPosib_start[l].first,j + southPosib_start[l].second };
						/*cellDetails[i][j].g = 0;
						cellDetails[i][j].h = 0;
						cellDetails[i][j].f++;*/
						cellDetails[i][j].parent = { i, j };
							
					}
					//std::cout << "south" << std::endl;
				}
				if (dir == direction::north)
				{
					neighbour = { i + northPosib[l].first,j + northPosib[l].second };
				}
				if (dir == direction::east)
				{
					neighbour = { i+eastPosib[l].first,j+eastPosib[l].second };
					//start variable
					if (i == src.first)
					{
						neighbour = { i + eastPosib_start[l].first,j + eastPosib_start[l].second };
						/*cellDetails[i][j].g = 0;
						cellDetails[i][j].h = 0;
						cellDetails[i][j].f++;*/
						cellDetails[i][j].parent = { i, j };
					}
					//std::cout << "east" << std::endl;
				}
				/*if (dir == direction::west)
				{
					neighbour = { i + westPosib[l].first,j + westPosib[l].second };
				}*/
				//std::cout << "neigh: " << neighbour.first << " " << neighbour.second << std::endl;
				//Pair neighbour(i + , j + add_y);
				// Only process this cell if this is a valid
				// one
				if (isValid(grid, neighbour)) {
					// If the destination cell is the same
					// as the current successor
					if (isDestination(neighbour, dest)) { // Set the Parent of // the destination cell

						cellDetails[neighbour.first][neighbour.second].parent = { i, j };
						//printf("The destination cell is found\n");
						tracePath(cellDetails, dest,res);
						return;
					}
					// If the successor is already on the
					// closed list or if it is blocked, then
					// ignore it. Else do the following
					else if (!closedList[neighbour.first][neighbour.second] && isUnBlocked(grid,neighbour))
					{
						double gNew, hNew, fNew;

						gNew = cellDetails[i][j].g + calculateGValue(grid, position, neighbour);
						if (dir == direction::south && (j == src.second))
						{
							//start variable
							gNew /= 2;
						}
						if (dir == direction::east && i == src.first)
						{
							//start variable
							gNew /= 2;
							
						}
						
						//hNew = calculateHValue(neighbour,
						//	dest);
						//hNew = 0;
						if (dir == direction::south)
						{
							//start variable
							if (j == src.second)
								hNew = 0;
							else
								hNew = abs(neighbour.first - dest.first) *0.01;
						}
						if (dir == direction::east)
						{
							//start variable
							if (i == src.first)
								hNew = 0;
							else
								hNew = abs(neighbour.second - dest.second)*0.01;
						}
						hNew = 0;
						//cout << "hnew: " << hNew << std::endl;
						//fNew = gNew + hNew*2;
						fNew = gNew + hNew;
						// If it isn’t on the open list, add
						// it to the open list. Make the
						// current square the parent of this
						// square. Record the f, g, and h
						// costs of the square cell
						//			 OR
						// If it is on the open list
						// already, check to see if this
						// path to that square is better,
						// using 'f' cost as the measure.
						if (cellDetails[neighbour.first]
							[neighbour.second].f == -1
							|| cellDetails[neighbour.first]
							[neighbour.second].f > fNew) 
						{
							openList.emplace(
								fNew, neighbour.first,
								neighbour.second);

							// Update the details of this
							// cell
							cellDetails[neighbour.first][neighbour.second].g = gNew;
							cellDetails[neighbour.first][neighbour.second].h = hNew;
							cellDetails[neighbour.first][neighbour.second].f = fNew;
							cellDetails[neighbour.first][neighbour.second].parent = { i, j };
							
		/*					if((dir == direction::south && j == src.second) || (dir == direction::east && i == src.first))
							{
								cellDetails[neighbour.first][neighbour.second].g = 0;
								cellDetails[neighbour.first][neighbour.second].h = 0;
								cellDetails[neighbour.first][neighbour.second].f = 0;
								cellDetails[neighbour.first][neighbour.second].parent = { neighbour.first, neighbour.second };
								openList.emplace(0.0, neighbour.first, neighbour.second);
							}*/
								
						}
					}
					
				}
			}
		}

		// When the destination cell is not found and the open
		// list is empty, then we conclude that we failed to
		// reach the destiantion cell. This may happen when the
		// there is no way to destination cell (due to
		// blockages)
		printf("Failed to find the Destination Cell\n");
	}
private:
	 


	// A structure to hold the necessary parameters
	struct cell {
		// Row and Column index of its parent
		Pair parent;
		// f = g + h
		double f, g, h;
		cell()
			: parent(-1, -1)
			, f(-1)
			, g(-1)
			, h(-1)
		{
		}
	};
	enum direction
	{
		south = 0,
		north = 1,
		west,
		east
	};


	// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
	
	bool isValid(const vector<vector<int>>& grid, const Pair& point)
	{ // Returns true if row number and column number is in
	// range
		if (ROW > 0 && COL > 0)
			return (point.first >= 0) && (point.first < ROW)
			&& (point.second >= 0)
			&& (point.second < COL);

		return false;
	}
	// A Utility Function to check whether the given cell is
	// blocked or not
	bool isUnBlocked(const vector<vector<int>>& grid,
		const Pair& point)
	{
		// Returns true if the cell is not blocked else false
		return isValid(grid, point)
			&& grid[point.first][point.second] != UINT8_MAX;
	}

	// A Utility Function to check whether destination cell has
	// been reached or not
	bool isDestination(const Pair& position, const Pair& dest)
	{
		return position == dest;
	}
	// A Utility Function to calculate the 'g' heuristics.
	
	double calculateGValue(const vector<vector<int>>& grid, const Pair& position, const Pair& neighbor)
	{
		double distance;
		if (position.first != neighbor.first && position.second != neighbor.second)
			distance = s_diagonal;
		else
			distance = 1;
		//distance = 1;
		return grid[neighbor.first][neighbor.second] * distance;
		
	}
	// A Utility Function to calculate the 'h' heuristics.
	double calculateHValue(const Pair& src, const Pair& dest)
	{
		// h is estimated with the two points distance formula
		/*return sqrt(pow((src.first - dest.first), 2.0)
			+ pow((src.second - dest.second), 2.0));*/
		return abs(src.first - dest.first) + abs(src.second - dest.second);
	}

	// A Utility Function to trace the path from the source to
	// destination
	
	void tracePath( const vector<vector<cell>>& cellDetails, const Pair& dest, vector<std::pair<int, int>>& res)
	{
	//	printf("\nThe Path is ");

		stack<Pair> Path;

		int row = dest.first;
		int col = dest.second;
		//Pair next_node = cellDetails[row][col].parent;
		Pair next_node = dest;
		
		do {
			//std::cout << row << " " << col << std::endl;
			Path.push(next_node);
			res.push_back(next_node);
			next_node = cellDetails[row][col].parent;
			row = next_node.first;
			col = next_node.second;
		} while (cellDetails[row][col].parent != next_node);
		Path.emplace(row, col);
		res.push_back(next_node);
		
		std::reverse(res.begin(), res.end());
	/*	vector<vector<int>> grid;
		grid.resize(ROW);
		for (int r = 0; r < ROW; r++)
			grid[r].resize(COL);*/
		//for (int o = 0; o < res.size(); o++)
			//std::cout << res[o].first <<" "<< res[o].second << std::endl;
		//std::cout << std::endl;
		//while (!Path.empty()) {
		//	Pair p = Path.top();
		//	std::cout << p.first << " " << p.second << std::endl;
		//	Path.pop();
		//	//grid[p.first][p.second] = 3;
		//	//printf("-> (%d,%d) ", p.first, p.second);
		//}
		//for (size_t i = 0; i < 60; i++)
		//{
		//	for (size_t j = 0; j < 60; j++)
		//	{
		//		//std::cout << grid[i][j] << " ";
		//	}
		//	//std::cout << std::endl;
		//}
		//std::cout << std::endl;
	}
	
	const double s_diagonal = sqrt(2);
	size_t ROW, COL;


};


#endif // !PATHFINDER_H
