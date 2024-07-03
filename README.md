# League of Legends Diamond Ranked Games (10 min) Dataset
https://www.kaggle.com/datasets/bobbyscience/league-of-legends-diamond-ranked-games-10-min/
## Overview
This project involves classifying League of Legends ranked games based on the first 10 minutes of gameplay using a neural network implemented from scratch in C++. The neural network aims to predict the outcome of high-ELO solo queue games.
## Achieved Accuracy
- The neural network achieved a mean accuracy of 72% in predicting the outcomes of the games using k-fold cross-validation.

  ![image](https://github.com/zainmo11/Lol-win-prediction/assets/89034348/6e769a85-8b9c-4617-b50d-a8111c502a9d)
## About the Dataset
### Context
League of Legends is a MOBA (Multiplayer Online Battle Arena) where two teams (Blue and Red) compete to destroy each other's Nexus. The map features three lanes, a jungle, and five roles. The primary goal is to destroy the enemy Nexus to win the game.

### Content
- Data: Contains statistics from approximately 10,000 ranked solo queue games after the first 10 minutes.
- Features: 38 features in total (19 per team) including kills, deaths, gold, experience, and level.
- Target: The blueWins column indicates if the Blue team won (1) or lost (0).
### Glossary
- Warding Totem: Reveals areas on the map.
- Minions: NPCs providing gold when killed.
- Jungle Minions: NPCs giving gold and buffs when killed.
- Elite Monsters: High HP/damage monsters providing significant bonuses when killed.
- Dragons: Elite monsters giving team-wide bonuses, with major advantages for the 4th and 5th dragons.
- Herald: Elite monster aiding in pushing lanes and destroying structures.
- Towers: Structures that need to be destroyed to reach the enemy Nexus.
- Level: Champion level, starting at 1 and maxing at 18.
## Architecture
### Neuron
The Neuron class represents a single neuron in the neural network. It encapsulates the following functionalities:

#### Data Members:

- m_nWeights: The number of weights associated with the neuron, including the bias weight.
- m_weights: A vector storing the weights of the neuron.
- m_activation: The weighted sum of inputs plus bias.
- m_output: The result of applying the activation function to the activation value.
- m_delta: The error term used for updating weights during backpropagation.
#### Methods:

- activate(const std::vector<float>& inputs): Computes the activation value by summing the weighted inputs plus bias.
- transfer(): Applies the sigmoid activation function to the activation value.
- transfer_derivative(): Computes the derivative of the sigmoid function, used for backpropagation.
- initWeights(int n_weights): Initializes weights with random values.
### Layer
The Layer class represents a layer of neurons in the neural network. It manages a collection of Neuron objects.

#### Data Members:

- m_neurons: A vector of Neuron objects in the layer.
#### Methods:

- initNeurons(int n_neurons, int n_weights): Initializes the neurons in the layer with the specified number of neurons and weights.
- get_neurons(): Provides access to the neurons in the layer.
### Network
The Network class represents the entire neural network and manages the layers, training, and prediction processes.

#### Data Members:

- m_nLayers: The number of layers in the network.
- m_layers: A vector of Layer objects representing the network's layers.
#### Methods:

- initialize_network(int n_inputs, int n_hidden, int n_outputs): Sets up the network with layers based on input, hidden, and output sizes.
- add_layer(int n_neurons, int n_weights): Adds a new layer to the network.
- forward_propagate(const std::vector<float>& inputs): Computes the network's output given a set of inputs.
- backward_propagate_error(const std::vector<float>& expected): Calculates error gradients and updates the delta values for neurons.
- update_weights(const std::vector<float>& inputs, float l_rate): Adjusts the weights based on the computed deltas and learning rate.
- train(const std::vector<std::vector<float>>& training_data, float l_rate, size_t n_epoch, size_t n_outputs): Trains the network using the provided training data.
- predict(const std::vector<float>& input): Predicts the output class for a given input.
- display_human() const: Displays the network's structure and neuron details in a human-readable format.


## Code Usage
### Main Function
The following code snippet demonstrates how to use the dataset with the provided neural network functions:

``` cpp
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <fstream>
#include <regex>
#include <algorithm>
#include <set>

std::vector<std::vector<float>> load_csv_data(std::string filename);
float accuracy_metric(std::vector<int> expect, std::vector<int> predict);
std::vector<float> evaluate_network(std::vector<std::vector<float>> dataset, int n_folds, float l_rate, int n_epoch, int n_hidden);

int main() {
    std::vector<std::vector<float>> csv_data = load_csv_data("../data.csv");

    // Normalize the last column for one-hot encoding
    std::map<int, int> lookup = {};
    int index = 0;
    for (auto& vec : csv_data) {
        auto ret = lookup.insert({static_cast<int>(vec.back()), index});
        vec.back() = static_cast<float>(ret.first->second);
        if (ret.second) {
            index++;
        }
    }

    int n_folds = 5;        // Number of folds for cross-validation
    float l_rate = 0.3f;    // Learning rate for the network
    int n_epoch = 500;      // Number of epochs for training
    int n_hidden = 5;       // Number of hidden neurons in the network

    // Evaluate the neural network
    std::vector<float> scores = evaluate_network(csv_data, n_folds, l_rate, n_epoch, n_hidden);

    // Calculate mean accuracy across folds
    float mean = std::accumulate(scores.begin(), scores.end(), 0.0f) / scores.size();
    std::cout << "Mean accuracy: " << mean << "%" << std::endl;

    return 0;
}
```
## Function Descriptions
- evaluate_network: Performs k-fold cross-validation on the dataset using a neural network. It returns accuracy scores for each fold.

```cpp
std::vector<float> evaluate_network(
    std::vector<std::vector<float>> dataset,
    int n_folds,
    float l_rate,
    int n_epoch,
    int n_hidden);
```
- accuracy_metric: Computes the accuracy of predictions compared to the expected results.

```cpp
float accuracy_metric(
    std::vector<int> expect,
    std::vector<int> predict);
```
- load_csv_data: Loads and normalizes data from a CSV file. The last column is left unchanged for target values.

```cpp
std::vector<std::vector<float>> load_csv_data(std::string filename);
```
## License
This code and dataset are provided under the MIT License. See the LICENSE file for details.
