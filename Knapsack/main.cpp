#include<bits/stdc++.h>
#include <algorithm>
#include <math.h>
#include <utility>
#include <stdio.h>


using namespace std;

const int POPULATION_SIZE = 4 ;
const int N_GENERATIONS = 100;
const double P_c = 0.7; //crossover rate
const double P_m =0.1; //mutation rate


//fitness = Summation from 1 to n of (--> chromosome * --> value) <= knapsack
vector<int> fitness (vector<int>& weight,vector<int>& value, vector<vector<int>>& population, int knapsack) {
    vector<int> fitness(POPULATION_SIZE);

    for (int i = 0; i < population.size(); ++i) {
        int sum1 = 0; //value
        int sum2 = 0; //weight

        for (int j = 0; j < population[i].size(); ++j) {
            sum1 += population[i][j] * value[j];
            sum2 += population[i][j] * weight[j];
        }

        if (sum2 <= knapsack) {
            fitness[i] = sum1;
        } else {
            fitness[i] = 0;
        }
    }
    return fitness;
}

vector<vector<int>> selection (vector<int>& fitness_v, int num_parents,vector<vector<int>>& population) {
    vector<int> tmp_f = fitness_v;
    vector<vector<int>> parents;

    for (int i = 0; i < num_parents; ++i) {
        auto maxIt = max_element(tmp_f.begin(), tmp_f.end());
        int maxIdx = distance(tmp_f.begin(), maxIt);

        parents.push_back(population[maxIdx]);
        *maxIt = -1e9;
    }
    return parents;
}

vector<vector<int>> crossover(vector<vector<int>>& parents, int n_offsprings) {
    vector<vector<int>> offsprings(n_offsprings, vector<int>(parents[0].size()));
    int r1 = parents[0].size() / 2;

    int i = 0;
    while (i < n_offsprings) {
        int parent1_i = i % parents.size();
        int parent2_i = (i + 1) % parents.size();
        double r2 = (double)(rand()) / RAND_MAX;

        if (r2 > P_c)
            continue;

        // Perform single point crossover
        for (int j = 0; j < r1; ++j)
            offsprings[i][j] = parents[parent1_i][j];

        for (int j = r1; j < parents[0].size(); ++j)
            offsprings[i][j] = parents[parent2_i][j];

        i++;
    }

    return offsprings;
}

void mutation(vector<vector<int>>& offsprings) {
    for (size_t i = 0; i < offsprings.size(); ++i) {
        for (size_t j = 0; j < offsprings[0].size(); ++j) {
            double random_value = (double)(rand()) / RAND_MAX;
            if (random_value < P_m) {
                // flip the bit 0->1 or 1->0
                offsprings[i][j] = 1 - offsprings[i][j];
            }
        }
    }
}


int main() {


    freopen("knapsack_input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int t; cin>>t;
    for(int i=1;i<=t;i++){
        cout << "Test Case " << i << endl;

        int knapsack,n_genes ;
        cin>>knapsack;
        cin>> n_genes;
        vector<int> weight(n_genes), value(n_genes);
        for(int i=0;i<n_genes;i++){
            cin>>weight[i]>>value[i] ;
        }


        vector<vector<int>>initial_population
                (POPULATION_SIZE, vector<int>(n_genes));
        for (int i = 0; i < POPULATION_SIZE; i++) {
            for (int j = 0; j < n_genes; j++) {
                initial_population[i][j] = rand() % 2;
            }
        }

        //select initial population
        cout << "Initial population: " << endl;
        for (int i = 0; i < POPULATION_SIZE; ++i)
        {
            for (int j = 0; j < n_genes; ++j)
            {
                cout << initial_population[i][j] << " ";
            }
            cout<<"c"<<i<<endl;
        }


//        vector<int> fitness_v;
//       fitness_v = fitness(weight, value, initial_population, knapsack_size);
//            for(auto it :fitness_v) {
//                cout<<it <<endl;
//            }


        vector<vector<int>> population = initial_population;

        for (int i = 0; i < N_GENERATIONS; ++i) {
            vector<int> fitness_v = fitness(weight, value, population, knapsack);
            vector<vector<int>> parents = selection(fitness_v, POPULATION_SIZE / 2, population);
            vector<vector<int>> offsprings = crossover(parents, POPULATION_SIZE - parents.size());
            mutation(offsprings);

            // Replace the old generation with new one
            population.clear();
            population.insert(population.end(), parents.begin(), parents.end());
            population.insert(population.end(), offsprings.begin(), offsprings.end());
        }

        // Find the best solution
        vector<int> last_generation_fitness = fitness(weight, value, population, knapsack);
        auto max_fitness = max_element(last_generation_fitness.begin(), last_generation_fitness.end());
        int idx = distance(last_generation_fitness.begin(), max_fitness);
        vector<int> chromosome = population[idx];


        // Print the selected chromosome
        cout << "Selected Chromosome:" << endl;
        for (size_t i = 0; i < chromosome.size(); ++i) {
            cout << chromosome[i] << " ";
        }
        cout << endl;

        // Calculate and print the number of selected items
        int numSelectedItems = 0;
        for (int i = 0; i < chromosome.size(); ++i) {
            if (chromosome[i] != 0) {
                numSelectedItems++;
            }
        }
        cout << "Number of Selected Items: " << numSelectedItems << endl;

        // Calculate and print the total value and total weight
        int totalValue = 0;
        int totalWeight = 0;
        for (int i = 0; i < chromosome.size(); ++i) {
            if (chromosome[i] != 0) {
                totalValue += value[i];
                totalWeight += weight[i];
            }
        }
        cout << "Total Value: " << totalValue << endl;
        cout << "Total Weight: " << totalWeight << endl;

        // Print the weight and value of each selected item
        cout << "Selected Items:" << endl;
        for (int i = 0; i < chromosome.size(); ++i) {
            if (chromosome[i] != 0) {
                cout << "Item " << i+1 << ": Weight = " << weight[i] << ", Value = " << value[i] << endl;

            }

        }
        cout<<endl;
        cout<<"--------------------------------------"<<endl;
        cout<<endl;



    }


    return 0;

}

