#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const int POPULATION_SIZE = 500;
const double Pm = 0.01;
const int MAX_GENERATIONS = 200;
const int Pc = 0.7;

struct Individual {
    vector<double> coefficients;
    double fitness;
};


double fitness(const vector<double>& coefficients, double x) {
    double result = 0.0;
    int D = coefficients.size() - 1;

    for (int i = D; i >= 0; --i)
        result += coefficients[i] * pow(x, D - i);

    return result;
}

//mean square error
double MSE(const vector<Individual>& population, const vector<pair<double, double>>& points) {
    double sum = 0.0;
    int N = points.size();
    for (const auto& individual : population) {
        for (const auto& point : points) {
            double predicted = fitness(individual.coefficients, point.first);
            double  tmp = predicted - point.second;
            sum += tmp * tmp;
        }
    }
    return sum / N ;
}

//Random population
vector<Individual> initializePopulation(int degree) {
    vector<Individual> population(POPULATION_SIZE);

    for (auto& individual : population) {
        individual.coefficients.resize(degree + 1);
        for (int i = 0; i <= degree; ++i) {
            individual.coefficients[i] = (rand() % 2000 - 1000) / 100.0; // Random values between -10 and 10
        }
    }

    return population;
}

// Tournament selection
Individual tournamentSelection(const vector<Individual>& population) {
    const int TOURNAMENT_SIZE =10;
    int randomIndex;
    Individual bestIndividual;

    for (int i = 0; i < TOURNAMENT_SIZE; ++i) {
        randomIndex = rand() % POPULATION_SIZE;
        if (i == 0 || population[randomIndex].fitness < bestIndividual.fitness) {
            bestIndividual = population[randomIndex];
        }
    }

    return bestIndividual;
}

// 2-point crossover
Individual crossover(const Individual& parent1, const Individual& parent2) {

    int crossoverPoint1 = rand() % parent1.coefficients.size();
    int crossoverPoint2 = rand() % parent1.coefficients.size();

    Individual child;
    child.coefficients.resize(parent1.coefficients.size());

    for (int i = 0; i < crossoverPoint1; ++i) {
        child.coefficients[i] = parent1.coefficients[i];
    }

    for (int i = crossoverPoint1; i < crossoverPoint2; ++i) {
        child.coefficients[i] = parent2.coefficients[i];
    }

    for (int i = crossoverPoint2; i < child.coefficients.size(); ++i) {
        child.coefficients[i] = parent1.coefficients[i];
    }

    return child;
}

// Non-uniform mutation
void mutate(Individual& individual) {
    for (auto& coefficient : individual.coefficients) {
        if ((rand() % 100) / 100.0 < Pm) {
            double Lx = coefficient - (-10); // delta lower = -10
            double Ux = 10 - coefficient;
            double r1 = (double) (rand()) / ((double)(RAND_MAX ));
            double y;
            if (r1 <= 0.5) {
                y = Lx;
            } else {
                y = Ux;
            }
            double r = (double) (rand()) / ((double)(RAND_MAX ));
            double delta = y * (1 - pow(r, 1 - static_cast<double>(1) / MAX_GENERATIONS));
            if (y == Lx){
                coefficient -= delta;
            } else {
                coefficient += delta;
            }

        }
    }

}

int main() {
    //input and output to files
    freopen("input.txt", "r", stdin);
    freopen("output2.txt", "w", stdout);

   //srand(static_cast<unsigned int>(time(0)));

    int numDatasets; cin >> numDatasets;

    for (int datasetIndex = 0; datasetIndex < numDatasets; ++datasetIndex) {
        int numPoints, degree;
        cin >> numPoints >> degree;

        vector<pair<double, double>> points(numPoints);
        for (int i = 0; i < numPoints; ++i) {
            cin >> points[i].first >> points[i].second;
        }

        vector<Individual> population = initializePopulation(degree);

        for (int generation = 0; generation < MAX_GENERATIONS; ++generation) {
            // Evaluate fitness
            for (auto& individual : population) individual.fitness = MSE({individual}, points);

            sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
                return a.fitness < b.fitness;});

            // Elitism: Keep the best individual
            Individual bestIndividual = population[0];

            // selection, crossover, and mutation
            vector<Individual> newPopulation(POPULATION_SIZE);
            newPopulation[0] = bestIndividual;

            for (int i = 1; i < POPULATION_SIZE; ++i) {
                Individual parent1 = tournamentSelection(population);
                Individual parent2 = tournamentSelection(population);
                Individual child = crossover(parent1, parent2);
                mutate(child);
                newPopulation[i] = child;
            }

            population = newPopulation;
        }

        cout << "Dataset " << datasetIndex+1 << endl;
         cout<< "Coefficients: ";
        for (const auto& coefficient : population[0].coefficients) {
            cout << coefficient << " ";
        }
        cout<<endl;
        cout << "MSE: " << population[0].fitness << endl;
        cout<<"----------------------------------------------" <<endl;
    }



    return 0;
}
