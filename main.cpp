#include <iostream>
#include <vector>
#include <iomanip>
#include "centrality_engine.h"

void print_centrality(const std::string &name, const std::vector<double> &scores)
{
    std::cout << "\n--- " << name << " ---" << std::endl;
    for (size_t i = 0; i < scores.size(); ++i)
    {
        if (scores[i] > 0 || name == "Betweenness Centrality" || name == "Degree Centrality")
        {
            std::cout << "Node " << std::setw(3) << i << ": " << std::fixed << std::setprecision(6) << scores[i] << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <graph_file.txt>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    try
    {
        std::cout << "Loading graph from " << filename << "..." << std::endl;
        Graph g = GraphLoader::loadFromEdgeList(filename);
        std::cout << "Graph loaded successfully with " << g.V << " vertices." << std::endl;

        auto degree_scores = Centrality::degree(g);
        print_centrality("Degree Centrality", degree_scores);

        auto closeness_scores = Centrality::closeness(g);
        print_centrality("Closeness Centrality", closeness_scores);

        auto betweenness_scores = Centrality::betweenness(g);
        print_centrality("Betweenness Centrality", betweenness_scores);

        auto eigenvector_scores = Centrality::eigenvector(g);
        print_centrality("Eigenvector Centrality", eigenvector_scores);

        auto pagerank_scores = Centrality::pageRank(g);
        print_centrality("PageRank", pagerank_scores);

        auto katz_scores = Centrality::katz(g);
        print_centrality("Katz Centrality", katz_scores);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}