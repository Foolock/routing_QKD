#include <iostream>
#include <vector>

using namespace std;

class Graph {
    int numVertices;
    vector<vector<int>> adjList;

public:
    Graph(int vertices) {
        numVertices = vertices;
        adjList.resize(numVertices);
    }

    void addEdge(int src, int dest) {
        adjList[src].push_back(dest);
        adjList[dest].push_back(src);
    }

    vector<vector<int>> getEdges() {
        vector<vector<int>> edges;
        for(int i=0; i<numVertices; i++) {
            for(int j=0; j<adjList[i].size(); j++) {
                int neighbor = adjList[i][j];
                if(i < neighbor) {
                    vector<int> edge;
                    edge.push_back(i);
                    edge.push_back(neighbor);
                    edges.push_back(edge);
                }
            }
        }
        return edges;
    }

    void printEdgesAsGrid() {
        vector<vector<int>> edges = getEdges();
        vector<vector<int>> grid(numVertices, vector<int>(numVertices, 0));
        for(int i=0; i<edges.size(); i++) {
            int src = edges[i][0];
            int dest = edges[i][1];
            grid[src][dest] = 1;
            grid[dest][src] = 1;
        }
        for(int i=0; i<numVertices; i++) {
            for(int j=0; j<numVertices; j++) {
                cout << grid[i][j] << " ";
            }
            cout << endl;
        }
    }
};

int main() {
    Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 4);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(2, 3);
    g.addEdge(3, 4);

    g.printEdgesAsGrid();

    return 0;
}
