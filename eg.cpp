
#include <vector>
#include "County.h"
#include <unordered_set>
using namespace std;

class CountyMap {
	private:
		///put stuff here
	public:
		//constructor
		//performance: O(n^2), n = number of counties
		CountyMap(vector<County> counties, vector<vector<bool>> adjacencyMatrix);

		CountyMap(const CountyMap& other);

		~CountyMap();

		//performance: O(1)
		County getCounty(int countyID);

		//performance: O(1)
		bool countyAdjacency(int countyID1, int countyID2);

		//performance: O(1)
		int getNumCounties();
		
		//performance: O(n), n = number of counties
		vector<int> getAdjacentCounties(int countyID);

		//performance: O(n), n = number of counties in the itinerary
		bool checkValidItinerary(vector<int> itinerary);
};

#endif
