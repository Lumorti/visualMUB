#include <SFML/Graphics.hpp>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>
#include <iomanip>
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

class Chain {
private:

	// Things needed regarding of visualisation
    int chainIndex = 0;
	int numVectors = 0;
    int dragging = -1;
	double radius = 0;
    std::vector<double> vectorLengths;
	double objective = 0.0;
	bool fixFirst = true;
	std::vector<double> angles;
	std::vector<std::pair<double,double>> positions;
	std::pair<double,double> startPosition;
	std::pair<int,int> vecIndices;
	std::vector<int> basisIndices;
	std::vector<std::pair<double, Chain*>> relation;
	std::vector<std::pair<double, int>> relationInds;
    double finalX = 0;
    double finalY = 0;
    double finalDist = 0;

	// Things only relevant for visualisation
	int lineThickness = 5;
    std::vector<sf::CircleShape> points;
    std::vector<sf::RectangleShape> lines;
	sf::CircleShape radiusCircle;
	sf::Text text;
	int pointSize = 5;

public:

	// Getter for the objective
	double getObjective() {
		return objective;
	}

    // Set whether the first angle should be fixed
    void setFirstAngleFixed(bool fixFirst_) {
        fixFirst = fixFirst_;
    }

    // Set the index
    void setIndex(int chainIndex_) {
        chainIndex = chainIndex_;
    }

    // Get the index
    int getIndex() {
        return chainIndex;
    }

	// Update the points and lines based on the angles
	void updatePointsAndLines() {

		// Update the points based on the angles
		std::pair<double,double> prevPos = startPosition;
		for (unsigned long int i = 0; i < points.size(); ++i) {
			positions[i] = {prevPos.first + vectorLengths[i] * sin(angles[i]), prevPos.second - vectorLengths[i] * cos(angles[i])};
			points[i].setPosition(positions[i].first, positions[i].second);
			prevPos = positions[i];
		}

		// Update the lines based on the points
		prevPos = startPosition;
		for (unsigned long int i = 0; i < lines.size(); ++i) {
			double angle = atan2(positions[i].second - prevPos.second, positions[i].first - prevPos.first);
			lines[i].setPosition(prevPos.first, prevPos.second);
			lines[i].setRotation(angle * 180.0 / M_PI);
            lines[i].setSize({float(vectorLengths[i]), float(lineThickness)});
			prevPos = positions[i];
		}

	}

	// Set how this chain is defined related to the others
	void setRelation(std::vector<std::pair<double, Chain*>> relation_, std::vector<std::pair<double, int>> relationInds_) {
		relation = relation_;
		relationInds = relationInds_;

		// Also change the colour to make it more obvious
		for (unsigned long int i = 0; i < points.size(); ++i) {
			points[i].setFillColor(sf::Color::Red);
		}

	}

    // Set the radius
    void setRadius(double radius_) {
        radius = radius_;
		radiusCircle.setRadius(radius);
		radiusCircle.setOrigin(radius, radius);
    }

    // Get the radius
    double getRadius() {
        return radius;
    }

    // Get the vector lengths
    std::vector<double> getVectorLengths() {
        return vectorLengths;
    }

    // Set the vector length
    void setVectorLengths(std::vector<double> vectorLengths_) {
        vectorLengths = vectorLengths_;
    }

    // Get a list of the chains that related to this
    std::vector<Chain*> getRelatedChains() {
        std::vector<Chain*> relatedChains;
        for (long unsigned int i = 0; i < relation.size(); ++i) {
            relatedChains.push_back(relation[i].second);
        }
        return relatedChains;
    }

    // Get the list of angles
    std::vector<double> getAngles() {
        return angles;
    }

    // Get the two vec indices
    std::pair<int, int> getVecIndices() {
        return vecIndices;
    }

    // Set the angles
    void setAngles(std::vector<double> angles_) {
        angles = angles_;
    }

	// Get the angle to the edge of the circle
	double getAngleToEdge() {
		double finalX = 0;
		double finalY = 0;
		for (unsigned long int i = 0; i < angles.size(); ++i) {
			finalX += vectorLengths[i] * sin(angles[i]);
			finalY -= vectorLengths[i] * cos(angles[i]);
		}
		return atan2(finalY, finalX);
	}

	// Snap the end point to the radius
	void snap() {

		// Repeat a few times for better accuracy
		for (long unsigned int i=0; i<10; i++) {

			// How close the last point is versus the circle radius
			double finalX = 0;
			double finalY = 0;
			for (unsigned long int i = 0; i < angles.size(); ++i) {
				finalX += vectorLengths[i] * sin(angles[i]);
				finalY -= vectorLengths[i] * cos(angles[i]);
			}
			double finalDist = radius - sqrt(finalX * finalX + finalY * finalY);
			double finalAngle = atan2(finalY, finalX);

			// If the final distance is less than the radius, snap it
			if (std::abs(finalDist) < 5.0) {

				// Move the last point this amount
				move({finalDist * cos(finalAngle), finalDist * sin(finalAngle)});

				// Update everything
				update();

			}

		}

	}

	// Check if the user is dragging any point
	bool isDragging() {
		return dragging != -1;
	}

	// Check if the user can modify this chain
	bool isFixed() {
		return relation.size() != 0;
	}

    // Would changing a certain angle of a certain chain affect the objective?
    double getGradient(int chainDiff, int angleDiff) {

        // If this chain is the one
        double coeff = 0.0;
        if (chainIndex == chainDiff) {
            coeff = 1.0;
        } 

        // Check if this chainDiff is part of the relation
        for (long unsigned int i=0; i<relationInds.size(); ++i) {
            if (relationInds[i].second == chainDiff) {
                coeff = relation[i].first;
                break;
            }
        }

        // Otherwise that chain and angle has no effect on this chain
        if (coeff == 0.0) {
            return 0.0;
        }

        // Otherwise calculate the derivative
        double P = finalDist - radius;
        double D = 1.0 / finalDist;
        double Z = coeff * 2.0 * vectorLengths[angleDiff] * (finalX * cos(angles[angleDiff]) + finalY * sin(angles[angleDiff]));
        double deriv = P * D * Z;
        return deriv;

    }

    // Would changing two angles of two chains affect the objective?
    double getGradient2(int chainDiff, int angleDiff, int chainDiff2, int angleDiff2) {

        // Whether the angle change causes a positive or negative rotation
        double coeff1 = 0.0;
        double coeff2 = 0.0;
        if (chainIndex == chainDiff) {
            coeff1 = 1.0;
        }
        if (chainIndex == chainDiff2) {
            coeff2 = 1.0;
        }
        for (long unsigned int i=0; i<relationInds.size(); ++i) {
            if (relationInds[i].second == chainDiff) {
                coeff1 = relationInds[i].first;
            }
            if (relationInds[i].second == chainDiff2) {
                coeff2 = relationInds[i].first;
            }
        }
        if (coeff1 == 0.0 && coeff2 == 0.0) {
            return 0.0;
        }

        double P = finalDist - radius;
        double D = 1.0 / finalDist;

        double L = vectorLengths[angleDiff];
        double Z1 = 2.0 * coeff1 * L * (finalX * cos(angles[angleDiff]) + finalY * sin(angles[angleDiff]));
        double Z2 = 2.0 * coeff2 * L * (finalX * cos(angles[angleDiff2]) + finalY * sin(angles[angleDiff2]));

        double dP = 0.5 * D * Z2;
        double dD = -0.5 * std::pow(D, 3) * Z2;
        double sinIfSame = 0;
        double cosIfSame = 0;
        if (chainDiff == chainDiff2 && angleDiff == angleDiff2) {
            sinIfSame = -sin(angles[angleDiff]);
            cosIfSame = cos(angles[angleDiff]);
        }
        double dZ1 = 2.0*L*L*coeff1*coeff2*cos(angles[angleDiff])*cos(angles[angleDiff2])
                    + 2*finalX*L*sinIfSame 
                    + 2.0*L*L*coeff1*coeff2*sin(angles[angleDiff])*sin(angles[angleDiff2])
                    + 2*finalY*L*cosIfSame; 

        double hessEl = dP*D*Z1 + P*(dD*Z1 + dZ1*D);

        //std::cout << "chain diff 1: " << chainDiff << std::endl;
        //std::cout << "angle diff 1: " << angleDiff << std::endl;
        //std::cout << "chain diff 2: " << chainDiff2 << std::endl;
        //std::cout << "angle diff 2: " << angleDiff2 << std::endl;
        //std::cout << "P: " << P << std::endl;
        //std::cout << "D: " << D << std::endl;
        //std::cout << "Z1: " << Z1 << std::endl;
        //std::cout << "Z2: " << Z2 << std::endl;
        //std::cout << "dP: " << dP << std::endl;
        //std::cout << "dD: " << dD << std::endl;
        //std::cout << "dZ1: " << dZ1 << std::endl;

        return hessEl;

    }

	// Update the objective
	void updateObjective() {

		// How close the last point is versus the circle radius
		finalX = 0;
		finalY = 0;
		for (unsigned long int i=0; i<angles.size(); ++i) {
			finalX += vectorLengths[i] * sin(angles[i]);
			finalY -= vectorLengths[i] * cos(angles[i]);
		}
		finalDist = sqrt(finalX * finalX + finalY * finalY);
		objective = std::pow(finalDist - radius, 2);

		// Set the label text
		text.setString(std::to_string(objective));

	}

	// Main constructor
    Chain(int numVectors_, std::pair<double,double> startPosition_, double vectorLength_, double circleRadius_, std::pair<int,int> vecIndices_, bool fixFirst_, std::vector<int> basisIndices_) {

		// Set the size of the points
		numVectors = numVectors_;
		vectorLengths = std::vector<double>(numVectors, vectorLength_);
		startPosition = startPosition_;
		vecIndices = vecIndices_;
		radius = circleRadius_;
		fixFirst = fixFirst_;
        basisIndices = basisIndices_;

		// For each vector
        for (int i = 0; i < numVectors; ++i) {

			// Create the point
            sf::CircleShape point(pointSize);
            point.setFillColor(sf::Color::Blue);
			point.setOrigin(pointSize, pointSize);
            points.push_back(point);

			// Create the angles and positions to be set later
			angles.push_back(0);
			positions.push_back(std::pair<double,double>(0,0));

			// Create the line
			sf::RectangleShape line(sf::Vector2f(vectorLengths[i], lineThickness));
			line.setFillColor(sf::Color::Black);
			line.setOrigin(0, lineThickness / 2.0);
			lines.push_back(line);

        }

		// Create the radius circle
		radiusCircle.setRadius(radius);
		radiusCircle.setOrigin(radius, radius);
		radiusCircle.setPosition(startPosition.first, startPosition.second);
		radiusCircle.setFillColor(sf::Color::Transparent);
		radiusCircle.setOutlineColor(sf::Color::Red);
		radiusCircle.setOutlineThickness(lineThickness);

		// Create the label
		text.setString(std::to_string(objective));
		text.setCharacterSize(20);
		text.setFillColor(sf::Color::Black);
		text.setPosition(startPosition.first, startPosition.second - radius - 20);

		// Set all the line positions
		update();

    }

    // Get the basis indices
    std::vector<int> getBasisIndices() {
        return basisIndices;
    }

	// Given a delta vector, move a given point
	void move(std::pair<double,double> delta, int toDrag=-1) {

		// If told to drag the last point
		if (toDrag == -1) {
			toDrag = points.size() - 1;
		}

		// The points after the dragged point
		for (unsigned long int i=toDrag; i<points.size(); ++i) {

			// Simple move with the point
			positions[i] = std::pair<double,double>(positions[i].first + delta.first, positions[i].second + delta.second);

		}

		// The points before the dragged point
		for (int i=toDrag-1; i>=0; --i) {

			// Get the angle between this point at the one after it
			double angle = atan2(positions[i + 1].second - positions[i].second, positions[i + 1].first - positions[i].first);
			
			// Move with that angle to the right distance
			positions[i] = std::pair<double,double>(positions[i + 1].first - vectorLengths[i + 1] * cos(angle), positions[i + 1].second - vectorLengths[i + 1] * sin(angle));

		}

		// Set all the angles
		std::pair<double,double> prevPos = startPosition;
		for (long unsigned int i=0; i<points.size(); ++i) {
			angles[i] = M_PI/2.0 + atan2(positions[i].second - prevPos.second, positions[i].first - prevPos.first);
			prevPos = positions[i];
		}

		// Fix the first angle
		if (fixFirst) {
			angles[0] = 0;
		}

	}

	// When an event happens, the chain should handle it
    bool handleEvent(const sf::Event& event, sf::RenderWindow& window) {

		// When the mouse button is pressed, the point being clicked should start being dragged
		if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {

			// Convert to coords
			sf::Vector2f mousePosition = window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y));

			// Check if the mouse is inside any of the points of this chain
			for (unsigned long int i = 0; i < points.size(); i++) {
				if (points[i].getGlobalBounds().contains(mousePosition)) {
					dragging = i;
					break;
				}
			}

			// If the mouse is inside any of the points, return true
			if (dragging != -1) {
				return true;
			}

		// When the mouse button is released, the point being dragged should stop
		} else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {

			// Only update if we were already dragging
			if (dragging != -1) {

				// Snap to the radius if nearby
				snap();

				// Stop dragging and update
				dragging = -1;
				return true;

			}

		// When the mouse is moved, the point being dragged should follow it
		} else if (event.type == sf::Event::MouseMoved) {

			// If there is a point being dragged
			if (dragging != -1) {

				// Convert to coords
				sf::Vector2f mousePosition = window.mapPixelToCoords(sf::Vector2i(event.mouseMove.x, event.mouseMove.y));

				// Move that point to the mouse, talking into account the size
				std::pair<double,double> delta = {mousePosition.x - positions[dragging].first, mousePosition.y - positions[dragging].second};
				move(delta, dragging);

				// Something changed so we should update everything
				return true;

			}

        }

		// Otherwise we don't need to update
		return false;

    }

	// Update the chain based on the others
    void updateFromRelations() {

        // If it's fixed, update the angles
        if (relation.size() != 0) {

            // For each element in the chain
            for (unsigned long int i = 0; i < points.size(); ++i) {

                // Get the sum of the angles from the other chains
                double angle = 0;
                for (unsigned long int j = 0; j < relation.size(); ++j) {
                    angle += relation[j].first * relation[j].second->angles[i];
                }

                // Set the angle to this point
                angles[i] = angle;

            }

        }

    }
    
	// Update the chain
	void update() {
		
        // If it's fixed, update the angles
        updateFromRelations();

		// Set all the point and line positions based on the angles
		updatePointsAndLines();

		// Update the objective
		updateObjective();

	}

	// Draw the points, lines, circle
    void draw(sf::RenderWindow& window, sf::Font& font) {
		window.draw(radiusCircle);
		for (const auto& line : lines) {
			window.draw(line);
		}
		for (unsigned long int i = 0; i < points.size(); ++i) {
            window.draw(points[i]);
        }
		text.setFont(font);
		window.draw(text);
    }

};

// Given a chain list, update them all and return the objective
double getObjective(std::vector<Chain>& chains) {

    // Update all the chains
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        chains[i].updateFromRelations();
        chains[i].updateObjective();
    }

    // The objective is the sum of the individual objectives
    double objective = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        objective += chains[i].getObjective();
    }

    // Return the objective
    return objective;

}

// Get the objective without updating anything
double getObjectiveNoUpdate(std::vector<Chain>& chains) {

    // The objective is the sum of the individual objectives
    double objective = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        objective += chains[i].getObjective();
    }

    // Return the objective
    return objective;

} 

// Given a list of chains, get the angles
std::vector<double> getAngles(std::vector<Chain>& chains) {

    // Get the angles from each chain
    std::vector<double> angles;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        if (chains[i].isFixed()) {
            continue;
        }
        std::vector<double> chainAngles = chains[i].getAngles();
        angles.insert(angles.end(), chainAngles.begin(), chainAngles.end());
    }

    // Return the angles
    return angles;

}

// Given a list of chains, get the angles
Eigen::VectorXd getAnglesEigen(std::vector<Chain>& chains) {

    // Get the angles from each chain
    std::vector<double> angles = getAngles(chains);
    Eigen::VectorXd anglesEigen(angles.size());
    for (unsigned long int i = 0; i < angles.size(); ++i) {
        anglesEigen(i) = angles[i];
    }

    // Return the angles
    return anglesEigen;

}

// Get the gradient of the objective at the current point
std::vector<double> getGradient(std::vector<Chain>& chains, bool fixFirst) {

    // Update all the chains
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        chains[i].updateFromRelations();
        chains[i].updateObjective();
    }

    // Differentiate by each angle in each chain
    std::vector<double> grad;
    int numAnglesPerChain = chains[0].getAngles().size();
    for (unsigned long int chainInd = 0; chainInd < chains.size(); chainInd++) {
        if (chains[chainInd].isFixed()) {
            continue;
        }
        for (int angleInd = 0; angleInd < numAnglesPerChain; angleInd++) {

            // Differentiate each chain by this angle
            double gradPerAngle = 0;
            if (!(fixFirst && angleInd == 0)) {
                for (unsigned long int i = 0; i < chains.size(); ++i) {
                    gradPerAngle += chains[i].getGradient(chainInd, angleInd);
                }
            }
            grad.push_back(gradPerAngle);

        }
    }

    // Return the objective
    return grad;

}

// Get the gradient of the objective at the current point
Eigen::VectorXd getGradientEigen(std::vector<Chain>& chains, bool fixFirst) {

    // Get the gradient as a vector
    std::vector<double> grad = getGradient(chains, fixFirst);

    // Convert to an Eigen vector
    Eigen::VectorXd gradEigen(grad.size());
    for (unsigned long int i = 0; i < grad.size(); ++i) {
        gradEigen(i) = grad[i];
    }

    // Return the objective
    return gradEigen;

}

// Get the Hessian of the objective at the current point
Eigen::MatrixXd getHessian(std::vector<Chain>& chains, bool fixFirst) {

    // Update all the chains
    int numAngles = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        chains[i].updateFromRelations();
        chains[i].updateObjective();
        if (!chains[i].isFixed()) {
            numAngles += chains[i].getAngles().size();
        }
    }

    // Differentiate by each angle in each chain
    int numAnglesPerChain = chains[0].getAngles().size();
    Eigen::MatrixXd hess(numAngles, numAngles);
    int ind1 = 0;
    for (unsigned long int chainInd = 0; chainInd < chains.size(); chainInd++) {
        if (chains[chainInd].isFixed()) {
            continue;
        }
        for (int angleInd = 0; angleInd < numAnglesPerChain; angleInd++) {

            // For the other chain and angle
            int ind2 = 0;
            for (unsigned long int chainInd2 = 0; chainInd2 < chains.size(); chainInd2++) {
                if (chains[chainInd2].isFixed()) {
                    continue;
                }
                for (int angleInd2 = 0; angleInd2 < numAnglesPerChain; angleInd2++) {

                    // Differentiate each chain by these two angles
                    double gradPerAngle = 0;
                    if (!(fixFirst && (angleInd == 0 || angleInd2 == 0))) {
                        for (unsigned long int i = 0; i < chains.size(); ++i) {
                            gradPerAngle += chains[i].getGradient2(chainInd, angleInd, chainInd2, angleInd2);
                        }
                    }
                    hess(ind1, ind2) = gradPerAngle;

                    ind2++;

                }
            }

            ind1++;

        }
    }

    // Return the Hessian
    return hess;

}

// Given a list of chains, set the angles
void setAngles(std::vector<Chain>& chains, std::vector<double> angles) {

    // Set the angles for each chain
    int index = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        if (chains[i].isFixed()) {
            continue;
        }
        std::vector<double> chainAngles(angles.begin() + index, angles.begin() + index + chains[i].getAngles().size());
        chains[i].setAngles(chainAngles);
        index += chains[i].getAngles().size();
    }

}

// Given a list of chains, set the angles
void setAngles(std::vector<Chain>& chains, Eigen::VectorXd angles) {

    // Set the angles for each chain
    int index = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        if (chains[i].isFixed()) {
            continue;
        }
        std::vector<double> chainAngles(angles.begin() + index, angles.begin() + index + chains[i].getAngles().size());
        chains[i].setAngles(chainAngles);
        index += chains[i].getAngles().size();
    }

}

// Data structure to pass along to the below function
struct ll_data_t {
    std::vector<Chain>* chains;
    bool fixFirst;
};

// The function used by the optim library
double optFunc(const Eigen::VectorXd& inputVals, Eigen::VectorXd* gradOut, void* optData) {

    // Extract the params needed
    ll_data_t* objfn_data = reinterpret_cast<ll_data_t*>(optData);
    std::vector<Chain>& chains = *(objfn_data->chains);
    bool fixFirst = objfn_data->fixFirst;

    // Set the angles
    std::vector<double> angles;
    for (long int i = 0; i < inputVals.size(); ++i) {
        angles.push_back(inputVals(i));
    }
    setAngles(chains, angles);

    // Update all the chains
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        chains[i].updateFromRelations();
        chains[i].updateObjective();
    }

    // The objective is the sum of the individual objectives
    double objective = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        objective += chains[i].getObjective();
    }

    // Differentiate by each angle in each chain
    if (gradOut) {
        int numAnglesPerChain = chains[0].getAngles().size();
        int ind = 0;
        for (unsigned long int chainInd = 0; chainInd < chains.size(); chainInd++) {
            if (chains[chainInd].isFixed()) {
                continue;
            }
            for (int angleInd = 0; angleInd < numAnglesPerChain; angleInd++) {

                // Differentiate each chain by this angle
                double gradPerAngle = 0;
                if (!(fixFirst && angleInd == 0)) {
                    for (unsigned long int i = 0; i < chains.size(); ++i) {
                        gradPerAngle += chains[i].getGradient(chainInd, angleInd);
                    }
                }
                (*gradOut)[ind] = gradPerAngle;
                ind++;

            }
        }
    }

    // Return the objective
    return objective;

}

// Entry point
int main(int argc, char* argv[]) {

	// Settings
	int d = -1;
	std::vector<int> N = {2, 2, 2};
	bool fixFirst = true;
    bool verbose = false;
    std::vector<std::string> modes;
	double startTemp = 1.0;
	int stepsFull = 100000;
	int stepsPartial = 1000;
    int maxIters = 1000;
    double lengthOG = 100;
    bool visual = true;

	// Parse the command line arguments
	std::string arg = "";
	for	(int i = 1; i < argc; ++i) {
		arg = argv[i];

		// Setting the basis sizes
		if (std::string(argv[i]) == "-N") {

			// Parse the comma seperated string
			N.clear();
			std::string NString = argv[i + 1];
			std::string delimiter = ",";
			size_t pos = 0;
			std::string token;
			while ((pos = NString.find(delimiter)) != std::string::npos) {
				token = NString.substr(0, pos);
				N.push_back(std::stoi(token));
				NString.erase(0, pos + delimiter.length());
			}
			N.push_back(std::stoi(NString));

		// If asked to display the help
		} else if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
			std::cout << "Usage: " << argv[0] << std::endl;
			std::cout << "General args: " << argv[0] << std::endl;
			std::cout << "  -N <basis sizes>   Comma seperated list of basis sizes, first is dimension" << std::endl; 
            std::cout << "  -1                 Don't fix the first angle" << std::endl;
            std::cout << "  -t <dbl>           Set the starting temperature" << std::endl;
            std::cout << "  -i <int>           Set the number of full anneal steps" << std::endl;
            std::cout << "  -I <int>           Set the number of partial anneal steps" << std::endl;
            std::cout << "  -p <int>           Set the number of times to anneal" << std::endl;
            std::cout << "  -d <int>           Set the dimension (changing only the radii)" << std::endl;
            std::cout << "  -V                 Verbose output" << std::endl;
            std::cout << "  -v                 Disable visualisation" << std::endl;
			std::cout << "Stackable args: " << argv[0] << std::endl;
            std::cout << "  --annealFull       Anneal all the chains at once" << std::endl;
            std::cout << "  --annealPartial    Anneal each chain one at a time" << std::endl;
            std::cout << "  --random           Randomise the angles" << std::endl;
            std::cout << "  --check            Check each chain in each direction" << std::endl;
            std::cout << "  --output           Output all angles in radians" << std::endl;
            std::cout << "  --outputDeg        Output all angles in degrees" << std::endl;
            std::cout << "  --outputTrue       Output the true angles" << std::endl;
            std::cout << "  --outputTrueDeg    Output the true angles in degrees" << std::endl;
            std::cout << "  --outputVectors    Output the true MU vectors" << std::endl;
            std::cout << "  --gradient         Perform gradient descent" << std::endl;
            std::cout << "  --optim            Use the optim library to minimize" << std::endl;
            std::cout << "  --shotgun          Run gradient descent from many points" << std::endl;
            std::cout << "  --optimshotgun     Run optim from many points" << std::endl;
            std::cout << "  --minima           Travel in the direction of decreasing local minima" << std::endl;
            std::cout << "  --decrease         Decrease the dimension by one, adiabatically" << std::endl;
            std::cout << "  --exact            Use known constructions to find an exact solution" << std::endl;
            std::cout << "  --stop             Stop, closing the window" << std::endl;
			return 0;

		// If told not to fix the first angle
		} else if (std::string(argv[i]) == "-1") {
			fixFirst = false;

        // If told to be verbose
        } else if (std::string(argv[i]) == "-V") {
            verbose = true;

        // If setting the starting temperature
        } else if (std::string(argv[i]) == "-t") {
            startTemp = std::stof(argv[i + 1]);

        // If told not to visualise
        } else if (std::string(argv[i]) == "-v") {
            visual = false;
            modes.push_back("stop");

        // If given a dimension
        } else if (std::string(argv[i]) == "-d") {
            d = std::stoi(argv[i + 1]);

        // If setting the anneal steps
        } else if (std::string(argv[i]) == "-i") {
            stepsFull = std::stoi(argv[i + 1]);

        // If setting the number iterations to anneal partial
        } else if (std::string(argv[i]) == "-I") {
            stepsPartial = std::stoi(argv[i + 1]);

        // If setting the number of times to anneal or shotgun
        } else if (std::string(argv[i]) == "-p") {
            maxIters = std::stoi(argv[i + 1]);

        // If given a word argument, add it the modes
        } else if (argv[i][0] == '-' && argv[i][1] == '-') { 
            modes.push_back(std::string(argv[i]).substr(2));

		}

	}

    // If there's a stop, move it to the end
    bool wasStop = false;
    for (unsigned long int i = 0; i < modes.size(); ++i) {
        if (modes[i] == "stop") {
            modes.erase(modes.begin() + i);
            wasStop = true;
            i--;
        }
    }
    if (wasStop) {
        modes.push_back("stop");
    }

    // If d hasn't been set, set it to the dimension of the first basis
    if (d == -1) {
        d = N[0];
    }

	// Create the main window
	sf::ContextSettings settings;
    settings.antialiasingLevel = 3.0;
	int windowWidth = 1280;
	int windowHeight = 720;
    if (!visual) {
        windowWidth = 1;
        windowHeight = 1;
    }
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Visual MUBs", sf::Style::Default, settings);

    if (visual) {
        window.setFramerateLimit(60);
    }

	// Load the font
	sf::Font font;
	if (!font.loadFromFile("/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf")) {
		std::cout << "Error loading font" << std::endl;
	}

	// Create each of the chains representing the MUBs
	double scaling = 100.0;
	double spacing = scaling * 3.0;
	double minX = 0;
	double minY = 0;
	double maxX = 0;
	double maxY = 0;
	std::vector<Chain> chains;
	for (unsigned long int i = 1; i < N.size(); ++i) {
		int NSoFarI = 0;
		for (unsigned long int m = 1; m < i; ++m) {
			NSoFarI += N[m];
		}
		for (unsigned long int j = i; j < N.size(); ++j) {
			int NSoFarJ = 0;
			for (unsigned long int m = 1; m < j; ++m) {
				NSoFarJ += N[m];
			}
			for (int k = 0; k < N[i]; ++k) {
				for (int l = 0; l < N[j]; ++l) {

					// The location of the chain
					int gridY = NSoFarI + k;
					int gridX = NSoFarJ + l;
					double currentX = gridX*spacing;
					double currentY = gridY*spacing;

					// Only the upper triangle
					if (gridX <= gridY) {
						continue;
					}

					// Update the min and max values
					minX = std::min(minX, currentX);
					minY = std::min(minY, currentY);
					maxX = std::max(maxX, currentX);
					maxY = std::max(maxY, currentY);
						
					// Orthogonality
					if (i == j) {
                        chains.push_back(Chain(N[0], {currentX, currentY}, sqrt(d)*scaling/double(d), 0.0, {gridX, gridY}, fixFirst, {int(i),int(k),int(j),int(l)}));

					// Mutually unbiasedness
					} else {
                        chains.push_back(Chain(N[0], {currentX, currentY}, sqrt(d)*scaling/double(d), scaling, {gridX, gridY}, fixFirst, {int(i),int(k),int(j),int(l)}));
					}

				}
			}
		}
	}

    // Assign an index to each chain
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        chains[i].setIndex(i);
    }

    // Create buttons to run the various commands
    std::vector<sf::RectangleShape> buttons;
    std::vector<sf::Text> buttonTexts;
    std::vector<std::string> buttonNames = {
                                    "output", 
                                    "outputDeg", 
                                    "outputTrue", 
                                    "outputTrueDeg",
                                    "outputVectors", 
                                    "random", 
                                    "check", 
                                    "check2", 
                                    "optim", 
                                    "optimshotgun", 
                                    "minima", 
                                    "decrease", 
                                    "exact", 
                                    "gradient", 
                                    "gradient2", 
                                    "shotgun", 
                                    "annealFull", 
                                    "annealPartial", 
                                    "fix/unfix first angle",
                                    "interrupt",
                                };
    for (unsigned long int i = 0; i < buttonNames.size(); ++i) {

        // Create the button
        sf::RectangleShape button(sf::Vector2f(200, 50));
        button.setFillColor(sf::Color::White);
        button.setOutlineColor(sf::Color::Black);
        button.setOutlineThickness(2.0);
        button.setPosition(minX - scaling*2 + 10, 10 + 60 * i);
        buttons.push_back(button);

        // Create the text for the button
        sf::Text buttonText;
        buttonText.setFont(font);
        buttonText.setString(buttonNames[i]);
        buttonText.setCharacterSize(20);
        buttonText.setFillColor(sf::Color::Black);
        buttonText.setPosition(minX - scaling*2 + 20, 20 + 60 * i);
        buttonTexts.push_back(buttonText);

    }

	// Create an empty eigen matrix
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(std::pow(chains.size(), 2), chains.size());
	int nextInd = 0;

	// Add the linear constraints
	for (unsigned long int i = 0; i < chains.size(); ++i) {
		std::pair<int,int> vecIndices1 = chains[i].getVecIndices();	
		for (unsigned long int j = i+1; j < chains.size(); ++j) {
			std::pair<int,int> vecIndices2 = chains[j].getVecIndices();

			// theta_abi + theta_cai = theta_cbi
			if (vecIndices1.first == vecIndices2.second) {

				// Find the other
				std::pair<int,int> lookingFor = std::make_pair(vecIndices2.first, vecIndices1.second);
				int otherInd = -1;
				for (unsigned long int k = 0; k < chains.size(); ++k) {
					if (chains[k].getVecIndices() == lookingFor) {
						otherInd = k;
						break;
					}
				}

				// Add to the matrix
				if (otherInd >= 0) {
                    if (chains.size() < 10) {
                        std::cout << "chain " << i << " + chain " << j << " = chain " << otherInd << std::endl;
                    }
					A(nextInd, i) = 1;
					A(nextInd, j) = 1;
					A(nextInd, otherInd) = -1;
					nextInd++;
				}

			}

		}
	}

	// Remove the zero rows
	A.conservativeResize(nextInd, Eigen::NoChange);

    // Print if it's not too big
    if (std::max(A.rows(), A.cols()) < 20 || verbose) {
        std::cout << "Original matrix:" << std::endl;
        std::cout << A << std::endl;
    }

    // Perform row reduction
    Eigen::MatrixXd reducedA = A;
	for (int i=0; i<std::min(reducedA.cols(), reducedA.rows()); ++i) {

        // Search for a non-zero element in this col
        int nonZeroInd = -1;
        Eigen::VectorXd nonZeroRow = Eigen::VectorXd::Zero(reducedA.cols());
        for (int j = i; j < reducedA.rows(); ++j) {
            if (std::abs(reducedA(j,i)) > 1e-10) {
                nonZeroInd = j;
                nonZeroRow = reducedA.row(j);
                break;
            }
        }
        std::cout << i << " " << nonZeroInd << std::endl;

        // If there is no non-zero element, skip this column
        if (nonZeroInd == -1) {
            continue;
        }

        // Make the corresponding row start with a one
        reducedA.row(i) += ((1.0-reducedA(i,i)) / nonZeroRow(i)) * nonZeroRow;

        // Zero all the other elements in this column
        for (int j = 0; j < reducedA.rows(); ++j) {
            if (j != i && std::abs(reducedA(j,i)) > 1e-10) {
                reducedA.row(j) -= (reducedA(j,i) / nonZeroRow(i)) * nonZeroRow;
            }
        }
        
    }

    // If the matrix isn't too big
    if (std::max(reducedA.rows(), reducedA.cols()) < 20 || verbose) {
        std::cout << "Reduced matrix:" << std::endl;
        std::cout << reducedA << std::endl;
    }

    // Calculate the maximum number of non-zero elements in any column
    int maxNonZero = 0;
    int minNonZero = 1000000;
    int avgNonZero = 0;
    for (int i = 0; i < reducedA.cols(); ++i) {
        int nonZero = 0;
        for (int j = 0; j < reducedA.rows(); ++j) {
            if (std::abs(reducedA(j,i)) > 1e-10) {
                nonZero++;
            }
        }
        maxNonZero = std::max(maxNonZero, nonZero);
        minNonZero = std::min(minNonZero, nonZero);
        avgNonZero += nonZero;
    }
    avgNonZero /= reducedA.cols();

	// For each row of the reduced matrix, get the first one and set the chains
	for (int i=0; i<std::min(reducedA.cols(), reducedA.rows()); ++i) {
        if (std::abs(reducedA(i,i)) > 1e-10) {
            std::vector<std::pair<double, Chain*>> terms;
            std::vector<std::pair<double, int>> termsInds;
            for (int k = i+1; k < reducedA.cols(); ++k) {
                if (std::abs(reducedA(i,k)) > 1e-10) {
                    terms.push_back(std::make_pair(-reducedA(i,k), &chains[k]));
                    termsInds.push_back(std::make_pair(-reducedA(i,k), k));
                }
            }
            chains[i].setRelation(terms, termsInds);
		}
	}

    // Count the number of chains with radius zero
    int numZero = 0;
    int numNonZero = 0;
    int numZeroAndFixed = 0;
    int numNonZeroAndFixed = 0;
    for (unsigned long int i = 0; i < chains.size(); ++i) {
        if (std::abs(chains[i].getRadius()) < 1e-10) {
            numZero++;
            if (chains[i].isFixed()) {
                numZeroAndFixed++;
            }
        } else {
            numNonZero++;
            if (chains[i].isFixed()) {
                numNonZeroAndFixed++;
            }
        }
    }
    int numFixed = numZeroAndFixed + numNonZeroAndFixed;
    int numFree = chains.size() - numFixed;

    // Output problem info
    std::cout << "dimension: " << d << std::endl;
    std::cout << "max non-zeros: " << maxNonZero << std::endl;
    std::cout << "min non-zeros: " << minNonZero << std::endl;
    std::cout << "avg non-zeros: " << avgNonZero << std::endl;
    std::cout << "we have " << chains.size() << " chains" << std::endl;
    std::cout << numFixed << " are fixed" << std::endl;
    std::cout << numFree << " are free" << std::endl;
    std::cout << numZero << " have radius zero" << std::endl;
    std::cout << "   " << numZeroAndFixed << " of these are fixed" << std::endl;
    std::cout << "   " << numZero - numZeroAndFixed << " of these are free" << std::endl;
    std::cout << chains.size() - numZero << " have radius non-zero" << std::endl;
    std::cout << "   " << numNonZeroAndFixed << " of these are fixed" << std::endl;
    std::cout << "   " << numNonZero - numNonZeroAndFixed << " of these are free" << std::endl;

    // Start with all random angles
    for (long unsigned int i=0; i<chains.size(); ++i) {
        std::vector<double> startingAngles = {0.0};
        for (int j=1; j<N[0]; ++j) {
            startingAngles.push_back(rand() / double(RAND_MAX) * 2.0 * M_PI);
        }
        chains[i].setAngles(startingAngles);
    }
    for (long unsigned int i=0; i<chains.size(); ++i) {
        chains[i].update();
    }

	// Global vars
	bool draggingBackground = false;
	sf::Vector2i lastMousePos;
	sf::View currentView = window.getView();
	sf::Vector2f lastViewPos;
	double zoomLevel = 1.0;
	double currentTemp = startTemp;
    double deltaTempFull = std::pow(1e-7 / startTemp, 1.0 / float(stepsFull));
    double deltaTempPartial = std::pow(1e-7 / startTemp, 1.0 / float(stepsPartial));
	double prevObjective = 100000000.0;
    int numDone = 0;
    double checkDelta = 3.00;
    sf::Clock clock;
    int fps = 0;
    int itersPerFrame = 1;
    int chainIndex = 0;
    std::string prevMode = "none";
    std::vector<double> allowedAngles;
    std::vector<double> bestMinimaAngles;
    double bestMinimaValue = 100000;
    std::vector<std::vector<double>> anglesToTest;
    int angleIndex = 0;
    double delMinima = 10.0;
    double alpha = 0.1;
    int numPerformed = 0;
    int numConverged = 0;
    double totalObjective = 0;
    double totalObjectiveLog = 0;

    // FPS and iters per draw counters
    sf::Text fpsCounter;
    fpsCounter.setFont(font);
    fpsCounter.setCharacterSize(20);
    fpsCounter.setFillColor(sf::Color::Black);
    fpsCounter.setPosition(-180, -100);
    sf::Text iterCounter;
    iterCounter.setFont(font);
    iterCounter.setCharacterSize(20);
    iterCounter.setFillColor(sf::Color::Black);
    iterCounter.setPosition(-180, -70);

	// Start the main loop
    while (window.isOpen() || !visual) {

        // Get the change in time, based on this do more iters per frame
        if (visual) {
            sf::Time elapsed = clock.getElapsedTime();
            clock.restart();
            fps = 1.0 / elapsed.asSeconds();
            if (fps > 10) {
                itersPerFrame += 5;
            } else if (itersPerFrame > 1) {
                itersPerFrame -= 1;
            }
            if (modes.size() == 0) {
                itersPerFrame = 0;
            } 
            fpsCounter.setString("FPS: " + std::to_string(fps));
            iterCounter.setString("iters per frame: " + std::to_string(itersPerFrame));
        } else {
            itersPerFrame = 100000000;
            fps = 0;
        }

		// Process events
        sf::Event event;
		bool somethingChanged = false;
        while (window.pollEvent(event)) {

			// Each chain should process events first
            for (auto& chain : chains) {
				bool val = chain.handleEvent(event, window);
				somethingChanged = somethingChanged || val;
			}

            // Check for button presses
			sf::Vector2f mousePosition = window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y));
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                for (unsigned long int i=0; i<buttons.size(); ++i) {
                    if (buttons[i].getGlobalBounds().contains(mousePosition)) {
                        std::cout << "pressed button: " << buttonNames[i] << std::endl;
                        if (buttonNames[i] == "fix/unfix first angle") {
                            fixFirst = !fixFirst;
                            for (auto& chain : chains) {
                                chain.setFirstAngleFixed(fixFirst);
                            }
                        } else if (buttonNames[i] == "interrupt") {
                            modes = {};
                        } else {
                            modes.push_back(buttonNames[i]);
                        }
                    }
                }
            }

			// When the window is closed, the program ends
            if (event.type == sf::Event::Closed) {
                window.close();

			// If we press the mouse
			} else if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {

				// Check if we are dragging any of the chains
				bool anyChainDragging = false;
				for (auto& chain : chains) {
					if (chain.isDragging()) { 
						anyChainDragging = true;
						break;
					}
				}

				// If not, we are dragging the background
				if (!anyChainDragging) {
					draggingBackground = true;
					lastMousePos = sf::Vector2i(event.mouseButton.x, event.mouseButton.y);
					lastViewPos = currentView.getCenter();
				}

			// If we release the mouse, we are no longer dragging the background
			} else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
				draggingBackground = false;

			// When we move the mouse while dragging the background
			} else if (event.type == sf::Event::MouseMoved && draggingBackground) {

				// Get the mouse movement since clicking
				sf::Vector2i currentMousePos = sf::Vector2i(event.mouseMove.x, event.mouseMove.y);
				sf::Vector2i delta = lastMousePos - currentMousePos;
				
				// Adjust the view
				currentView.setCenter(lastViewPos.x + zoomLevel*delta.x, lastViewPos.y + zoomLevel*delta.y);
				window.setView(currentView);

			// When zooming with the mouse wheel
			} else if (event.type == sf::Event::MouseWheelScrolled) {
				currentView.zoom(1.0 - event.mouseWheelScroll.delta*0.1);
				zoomLevel *= 1.0 - event.mouseWheelScroll.delta*0.1;
				window.setView(currentView);

			// When resized
			} else if (event.type == sf::Event::Resized) {
				currentView.setSize(event.size.width, event.size.height);
				window.setView(currentView);
            }
            
        }

        // Get the next mode
        std::string mode = "none";
        if (modes.size() > 0) {
            mode = modes[0];
        }
        if (mode != prevMode) {
            prevMode = mode;
            itersPerFrame = 1;
            currentTemp = startTemp;
            numDone = 0;
            checkDelta = 3.0;
            alpha = 0.1;
            numPerformed = 0;
            numConverged = 0;
            delMinima = 10.0;
            std::cout << std::endl;
        }

        // If told to stop
        if (mode == "stop") {
            return 0;

		// If told to anneal
        } else if (mode == "annealFull") {

            // Do a bunch of iters at once
            double overallObjective = 0.0;
            double averageObjective = 0.0;
            double maxObjective = 0.0;
            int numAccepts = 0;
            for (int iter = 0; iter < itersPerFrame; ++iter) {

                // For each chain
                std::vector<std::vector<double>> prevAngles(chains.size());
                for (unsigned long int i = 0; i < chains.size(); ++i) {

                    // Wiggle each angle a bit
                    std::vector<double> angles = chains[i].getAngles();
                    prevAngles[i] = angles;
                    if (chains[i].isFixed()) {
                        continue;
                    }
                    for (unsigned long int j = 0; j < angles.size(); ++j) {
                        angles[j] += currentTemp * ((rand() / (double)RAND_MAX)-0.5);
                    }
                    if (fixFirst) {
                        angles[0] = 0.0;
                    }
                    chains[i].setAngles(angles);

                }

                // See if the overall objective function is better
                overallObjective = 0.0;
                maxObjective = 0.0;
                averageObjective = 0.0;
                for (auto& chain : chains) {
                    chain.update();
                }
                for (auto& chain : chains) {
                    overallObjective += std::pow(chain.getObjective(), 2);
                    maxObjective = std::max(maxObjective, chain.getObjective());
                    averageObjective += chain.getObjective();
                }
                averageObjective /= chains.size();

                // Accept with Boltzmann probability
                double randVal = rand() / double(RAND_MAX);
                double relativeDiff = 100.0 * (prevObjective - overallObjective) / prevObjective;
                double prob = std::min(1.0, std::exp(relativeDiff / currentTemp));
                if (randVal < prob) {
                    numAccepts++;
                    prevObjective = overallObjective;
                } else {
                    for (unsigned long int i = 0; i < chains.size(); ++i) {
                        chains[i].setAngles(prevAngles[i]);
                    }
                }

                // Lower the temperature
                currentTemp *= deltaTempFull;
                if (currentTemp < 1e-6) {
                    break;
                }

            }

            // Per iteration output
            std::cout << "sqr=" << overallObjective << "  avg=" << averageObjective << "  max=" << maxObjective << "  tmp=" << currentTemp << "  acp=" << 100.0 * float(numAccepts) / float(itersPerFrame) << "       \r" << std::flush;

            // If we are done, remove this mode
            if (currentTemp < 1e-6) {
                currentTemp = startTemp;
                modes.erase(modes.begin());
            }

        // If told to anneal little sections at once 
        } else if (mode == "annealPartial") {

            // Do a bunch of iters
            std::vector<Chain*> relatedChains = chains[chainIndex].getRelatedChains();
            double overallObjective = 0.0;
            double maxObjective = 0.0;
            double averageObjective = 0.0;
            for (int iter=0; iter<itersPerFrame; ++iter) {

                // For each chain
                std::vector<std::vector<double>> prevAngles(relatedChains.size());
                for (unsigned long int i=0; i<relatedChains.size(); ++i) {

                    // Wiggle each angle a bit
                    std::vector<double> angles = relatedChains[i]->getAngles();
                    prevAngles[i] = angles;
                    for (unsigned long int j=0; j<angles.size(); ++j) {
                        angles[j] += currentTemp * ((rand() / (double)RAND_MAX)-0.5);
                    }
                    if (fixFirst) {
                        angles[0] = 0.0;
                    }
                    relatedChains[i]->setAngles(angles);

                }

                // See if the overall objective function is better
                overallObjective = 0.0;
                maxObjective = 0.0;
                averageObjective = 0.0;
                for (auto& chain : chains) {
                    chain.update();
                }
                for (auto& chain : chains) {
                    overallObjective += std::pow(chain.getObjective(), 2);
                    maxObjective = std::max(maxObjective, chain.getObjective());
                    averageObjective += chain.getObjective();
                }
                averageObjective /= chains.size();

                // Accept with Boltzmann probability
                double randVal = rand() / (double)RAND_MAX;
                double prob = std::exp((prevObjective - overallObjective) / currentTemp);
                if (randVal < prob) {
                    prevObjective = overallObjective;
                } else {
                    for (unsigned long int i=0; i<relatedChains.size(); ++i) {
                        relatedChains[i]->setAngles(prevAngles[i]);
                    }
                }

                // Lower the temperature
                currentTemp *= deltaTempPartial;

                // Stop when we get to a low enough temperature
                if (currentTemp < 1e-6) {
                    numDone++;
                    currentTemp = startTemp;
                    chainIndex = rand() % chains.size();
                    while (!chains[chainIndex].isFixed()) {
                        chainIndex = rand() % chains.size();
                    }
                    break;
                }
                
            }

            // Per iteration output
            std::cout << numDone << "  sqr=" << overallObjective << "  avg=" << averageObjective << "  max=" << maxObjective << "  tmp=" << currentTemp << "        \r" << std::flush;

            // If we're done
            if (numDone > maxIters) {
                numDone = 0;
                modes.erase(modes.begin());
            }

        // If told to randomise the angles
        } else if (mode == "random") {

            // For each chain
            for (auto& chain : chains) {

                // Set each angle to a random value
                std::vector<double> angles = chain.getAngles();
                for (unsigned long int j=0; j<angles.size(); ++j) {
                    angles[j] = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
                }
                if (fixFirst) {
                    angles[0] = 0.0;
                }
                chain.setAngles(angles);

            }

            // Update the chains
            for (auto& chain : chains) {
                chain.update();
            }

            // If we're done
            modes.erase(modes.begin());

        // If setting to exact solution
        } else if (mode == "exact") {

            std::cout << "finding exact solution" << std::endl;

            // Number of bases
            int n = N.size();

            // Create the X operator
            Eigen::MatrixXcd X = Eigen::MatrixXcd::Zero(d,d);
            for (int i=0; i<d-1; i++) {
                X(i+1,i) = 1;
            }
            X(0,d-1) = 1;

            // Create the Z operator
            Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(d,d);
            std::complex<double> omega = std::exp(2.0*M_PI*std::complex<double>(0,1)/double(d));
            for (int i=0; i<d; i++) {
                Z(i,i) = std::pow(omega, i);
            }

            // Create the various operators (X, Z, XZ, XZZ, etc.)
            std::vector<Eigen::MatrixXcd> ops = {Z, X};
            for (int i=0; i<n-2; i++) {
                Eigen::MatrixXcd newOp = Z;
                for (int j=0; j<i; j++) {
                    newOp = Z*newOp;
                }
                ops.push_back(X*newOp);
            }

            // Get the eigenbases of each operator
            std::vector<std::vector<Eigen::VectorXcd>> eigenbases;
            for (int i=0; i<ops.size(); i++) {
                std::cout << "Basis " << i << ":" << std::endl;
                Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(ops[i]);
                std::vector<Eigen::VectorXcd> newBasis;
                for (int j=0; j<d; j++) {
                    newBasis.push_back(es.eigenvectors().col(j));
                    std::cout << "Vector " << j << ": " << newBasis[j].transpose() << std::endl;
                    std::cout << "Eigenvalue: " << es.eigenvalues()[j] << std::endl;
                }
                std::cout << std::endl << std::endl;
                eigenbases.push_back(newBasis);
            }

            // For each vector in each basis
            double maxError = 0;
            for (int i=0; i<n; i++) {
                for (int k=0; k<N[i]; k++) {

                    // For each other vector in each basis
                    for (int j=i+1; j<n; j++) {
                        for (int l=0; l<N[j]; l++) {

                            // Calculate the error
                            maxError = std::max(maxError, std::abs(std::abs(eigenbases[i][k].dot(eigenbases[j][l]))-1.0/std::sqrt(double(d))));

                        }
                    }

                }
            }
            std::cout << "d=" << d << " n=" << n << " error=" << maxError << std::endl;

            // Put these values into the chains
            for (auto& chain : chains) {
                std::pair<int,int> vecIndices = chain.getVecIndices();
                std::vector<double> angles(d);
                for (int i=0; i<d; i++) {
                    std::vector<int> inds = chain.getBasisIndices();
                    double angle1 = std::arg(eigenbases[inds[0]][inds[1]][i]);
                    double angle2 = std::arg(eigenbases[inds[2]][inds[3]][i]);
                    angles[i] = angle1 - angle2;
                }
                chain.setAngles(angles);
            }

            // Update the chains
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].update();
                }
            }

            // Done
            modes.erase(modes.begin());


        // If decreasing the dimension TODO
        } else if (mode == "decrease") {

            for (int l=0; l<itersPerFrame; l++) {

                // Start from the current angles
                Eigen::VectorXd x = getAnglesEigen(chains);

                // Settings for optim
                optim::algo_settings_t settings;
                settings.print_level = verbose ? 1 : 0;
                settings.iter_max = 100;

                // The data object passed to the optim function
                ll_data_t opt_data;
                opt_data.chains = &chains;
                opt_data.fixFirst = fixFirst;

                // Run the optimisation
                bool success = optim::lbfgs(x, optFunc, &opt_data, settings);

                // If any of the values are infinite, skip
                bool allGood = true;
                for (int i=0; i<x.size(); i++) {
                    if (std::abs(x(i)) > 1e10) {
                        allGood = false;
                        break;
                    }
                }
                if (!allGood) {
                    std::cout << "Some values are infinite, stopping" << std::endl;
                    break;
                }

                // Run the final values
                setAngles(chains, x);
                double finalObj = getObjective(chains);

                // Adjust the lengths
                double deltaL = 0.001;
                double ogLength = sqrt(d)*scaling/double(d);
                double refLength = 100;
                unsigned long int numDone = 0;
                double dPrime = d;
                for (auto& chain : chains) {
                    std::vector<double> lengths = chain.getVectorLengths();
                    for (unsigned long int i=0; i<lengths.size(); ++i) {
                        if (int(i) == d-1) {
                            lengths[i] -= deltaL;
                            if (lengths[i] < 0.0) {
                                lengths[i] = 0.0;
                                numDone++;
                            }
                            refLength = lengths[i];
                        }
                    }
                    chain.setVectorLengths(lengths);
                    double currentRad = chain.getRadius();
                    dPrime = double(d) + (refLength/ogLength) - 1.0;
                    if (std::abs(currentRad) > 1e-10) {
                        currentRad = dPrime * ogLength / sqrt(dPrime);
                        chain.setRadius(currentRad);
                    }
                }
                std::cout << "d = " << dPrime << " obj = " << finalObj << "  iters = " << itersPerFrame << "  \r" << std::flush;
            }
            std::cout << std::endl;

            // If all chains at the right length, stop
            if (numDone == chains.size()) {
                modes.erase(modes.begin());
            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // If using the check mode (checking each chain in each direction)
        } else if (mode == "check") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            double maxDiff = 0.0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                maxDiff = std::max(maxDiff, chains[i].getObjective());
            }
            std::cout << "sqr=" << bestObjective << "  max=" << maxDiff << "  del=" << checkDelta << "     \r" << std::flush;

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Keep track of each time we update the best objective
                int numChanges = 0;

                // For each non-fixed chain
                for (unsigned long int j=0; j<chains.size(); ++j) {
                    if (chains[j].isFixed()) {
                        continue;
                    }

                    // Check the effect of moving each direction
                    std::vector<double> angles = chains[j].getAngles();
                    for (unsigned long int k=0; k<angles.size(); ++k) {
                        if (k == 0 && fixFirst) {
                            continue;
                        }

                        // Save the old angle
                        double oldAngle = angles[k];

                        // Check plus checkDelta
                        angles[k] += checkDelta;
                        chains[j].setAngles(angles);
                        double objective = getObjective(chains);
                        if (objective >= bestObjective) {
                            angles[k] = oldAngle;
                        } else {
                            bestObjective = objective;
                            numChanges++;
                            continue;
                        }

                        // Check minus checkDelta
                        angles[k] -= checkDelta;
                        chains[j].setAngles(angles);
                        objective = getObjective(chains);
                        if (objective >= bestObjective) {
                            angles[k] = oldAngle;
                        } else {
                            bestObjective = objective;
                            numChanges++;
                            continue;
                        }

                    }
                    chains[j].setAngles(angles);

                }

                // If no changes, decrease checkDelta
                if (numChanges == 0) {
                    checkDelta *= 0.1;

                    // If it's small enough, stop
                    if (checkDelta < 1e-30) {
                        checkDelta = 3.0;
                        modes.erase(modes.begin());
                        break;
                    }
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // If using the check mode (checking each chain in each direction)
        } else if (mode == "check2") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            double maxDiff = 0.0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                maxDiff = std::max(maxDiff, chains[i].getObjective());
            }
            std::cout << "sqr=" << bestObjective << "  max=" << maxDiff << "  del=" << checkDelta << "     \r" << std::flush;

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Keep track of each time we update the best objective
                int numChanges = 0;

                // For each non-fixed chain
                for (unsigned long int j=0; j<chains.size(); ++j) {
                    if (chains[j].isFixed()) {
                        continue;
                    }

                    // Check the effect of moving each direction
                    std::vector<double> angles = chains[j].getAngles();
                    for (unsigned long int k=0; k<angles.size(); ++k) {
                        for (unsigned long int k2=0; k2<angles.size(); ++k2) {
                            if (k == 0 && fixFirst) {
                                continue;
                            }
                            if (k == k2) {
                                continue;
                            }

                            // Save the old angle
                            double oldAngle = angles[k];
                            double oldAngle2 = angles[k2];

                            // Check plus plus
                            angles[k] += checkDelta;
                            angles[k2] += checkDelta;
                            chains[j].setAngles(angles);
                            double objective = getObjective(chains);
                            if (objective >= bestObjective) {
                                angles[k] = oldAngle;
                                angles[k2] = oldAngle2;
                            } else {
                                bestObjective = objective;
                                numChanges++;
                                continue;
                            }

                            // Check plus minus
                            angles[k] += checkDelta;
                            angles[k2] -= checkDelta;
                            chains[j].setAngles(angles);
                            objective = getObjective(chains);
                            if (objective >= bestObjective) {
                                angles[k] = oldAngle;
                                angles[k2] = oldAngle2;
                            } else {
                                bestObjective = objective;
                                numChanges++;
                                continue;
                            }

                            // Check minus plus
                            angles[k] -= checkDelta;
                            angles[k2] += checkDelta;
                            chains[j].setAngles(angles);
                            objective = getObjective(chains);
                            if (objective >= bestObjective) {
                                angles[k] = oldAngle;
                                angles[k2] = oldAngle2;
                            } else {
                                bestObjective = objective;
                                numChanges++;
                                continue;
                            }

                            // Check minus minus
                            angles[k] -= checkDelta;
                            angles[k2] -= checkDelta;
                            chains[j].setAngles(angles);
                            objective = getObjective(chains);
                            if (objective >= bestObjective) {
                                angles[k] = oldAngle;
                                angles[k2] = oldAngle2;
                            } else {
                                bestObjective = objective;
                                numChanges++;
                                continue;
                            }

                        }
                    }
                    chains[j].setAngles(angles);

                }

                // If no changes, decrease checkDelta
                if (numChanges == 0) {
                    checkDelta *= 0.1;

                    // If it's small enough, stop
                    if (checkDelta < 1e-30) {
                        checkDelta = 3.0;
                        modes.erase(modes.begin());
                        break;
                    }
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Test many random starts followed by gradient descent
        } else if (mode == "shotgun") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            std::cout << std::setprecision(4);
            std::cout << "sqr=" << bestObjective << "  alp=" << alpha << "  con= " << numConverged << " / " << numPerformed << "  " << chains.size() << " " << numFree << " " << numFree*d << "     \r" << std::flush;

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Get the gradient
                std::vector<double> currentAngles = getAngles(chains);
                std::vector<double> gradient = getGradient(chains, fixFirst);
                double newObj = getObjectiveNoUpdate(chains);

                // Check if we need to decrease alpha
                if (std::isnan(newObj) || newObj > bestObjective || std::abs(newObj - bestObjective) < 1e-20) {
                    alpha *= 0.8;
                } else {
                    bestObjective = newObj;
                }

                // Update the angles based on this
                for (unsigned long int i=0; i<currentAngles.size(); ++i) {
                    currentAngles[i] -= gradient[i] * alpha;
                }
                setAngles(chains, currentAngles);

                // Stop when the gradient is nothing
                if (alpha < 1e-15) {

                    // Keep track of the successes versus attempts
                    if (bestObjective < 1e-7) {
                        numConverged++;
                    }
                    numPerformed++;
                    if (!std::isnan(bestObjective)) {
                        totalObjective += bestObjective;
                        totalObjectiveLog += std::log10(bestObjective);
                    }
                    std::cout << "s=" << bestObjective << " a=" << alpha << " c= " << numConverged << " / " << numPerformed << "  " << chains.size() << " " << numFree << " " << numFree*d << "     \n" << std::flush;
                    if (numPerformed >= maxIters) {
                        modes.erase(modes.begin());
                        break;
                    }

                    // Otherwise start again with random angles
                    alpha = 0.1;
                    bestObjective = 1e10;
                    for (auto& chain : chains) {
                        std::vector<double> angles = chain.getAngles();
                        for (unsigned long int j=0; j<angles.size(); ++j) {
                            angles[j] = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
                        }
                        if (fixFirst) {
                            angles[0] = 0.0;
                        }
                        chain.setAngles(angles);
                    }

                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Perform gradient descent
        } else if (mode == "gradient") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            double maxDiff = 0.0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                maxDiff = std::max(maxDiff, chains[i].getObjective());
            }
            std::cout << "sqr=" << bestObjective << "  max=" << maxDiff << "  alpha=" << alpha << "     \r" << std::flush;

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Get the gradient
                std::vector<double> currentAngles = getAngles(chains);
                Eigen::VectorXd gradient = getGradientEigen(chains, fixFirst);
                double newObj = getObjectiveNoUpdate(chains);

                // Decrease alpha if we need to
                if (std::isnan(newObj) || newObj > bestObjective || gradient.norm() < 1e-20 || std::abs(newObj - bestObjective) < 1e-20) {
                    alpha *= 0.8;
                } else {
                    bestObjective = newObj;
                }

                // Update the angles based on this
                for (unsigned long int i=0; i<currentAngles.size(); ++i) {
                    currentAngles[i] -= gradient[i] * alpha;
                }
                setAngles(chains, currentAngles);

                // Stop when the gradient is nothing
                if (alpha < 1e-10) {
                    modes.erase(modes.begin());
                    break;
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Perform gradient descent with the Hessian
        } else if (mode == "gradient2") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            double maxDiff = 0.0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                maxDiff = std::max(maxDiff, chains[i].getObjective());
            }
            std::cout << "sqr=" << bestObjective << "  max=" << maxDiff << "  alpha=" << alpha << "     \r" << std::flush;

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Get the gradient and hessian
                std::vector<double> currentAngles = getAngles(chains);
                Eigen::MatrixXd hessian = getHessian(chains, fixFirst);
                Eigen::VectorXd gradient = getGradientEigen(chains, fixFirst);
                double newObj = getObjectiveNoUpdate(chains);

                // Check if it's PSD
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hessian);
                Eigen::VectorXd eigenvalues = es.eigenvalues();
                double minEigenvalue = eigenvalues.minCoeff();

                // Get the search direction
                Eigen::VectorXd direction = -gradient;
                if (minEigenvalue >= 0.0) {
                    direction = hessian.colPivHouseholderQr().solve(-gradient);
                }

                std::cout << "hessian=\n" << hessian << std::endl;
                std::cout << "min eigenvalue=" << minEigenvalue << std::endl;
                std::cout << "gradient=\n" << gradient << std::endl;
                std::cout << "direction=\n" << direction << std::endl;
                modes.erase(modes.begin());
                break;

                // Decrease alpha if we need to
                if (std::isnan(newObj) || newObj > bestObjective - 1e-10) {
                    alpha *= 0.8;
                } else {
                    bestObjective = newObj;
                }

                // Update the angles based on this
                for (unsigned long int i=0; i<currentAngles.size(); ++i) {
                    currentAngles[i] += direction[i] * alpha;
                }
                setAngles(chains, currentAngles);

                // Stop when the gradient is nothing
                if (alpha < 1e-10) {
                    modes.erase(modes.begin());
                    break;
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Use the optim library
        } else if (mode == "optim") {

            // Only do this once
            modes.erase(modes.begin());

            // Start from the current angles
            Eigen::VectorXd x = getAnglesEigen(chains);

            // Settings for optim
            optim::algo_settings_t settings;
            settings.print_level = 1;
            settings.iter_max = 10000000;

            // The data object passed to the optim function
            ll_data_t opt_data;
            opt_data.chains = &chains;
            opt_data.fixFirst = fixFirst;

            // Run the optimisation
            bool success = optim::lbfgs(x, optFunc, &opt_data, settings);
            std::cout << "success = " << success << std::endl;

            // Run the final values
            setAngles(chains, x);
            double finalObj = getObjective(chains);
            std::cout << "finalObj = " << finalObj << std::endl;

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Use the optim library
        } else if (mode == "optimshotgun") {

            // Calculate the objective
            double bestObjective = getObjective(chains);
            std::cout << std::setprecision(4);

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Start from the current angles
                Eigen::VectorXd x = getAnglesEigen(chains);

                // Settings for optim
                optim::algo_settings_t settings;
                settings.print_level = 0;
                settings.iter_max = 1000;

                // The data object passed to the optim function
                ll_data_t opt_data;
                opt_data.chains = &chains;
                opt_data.fixFirst = fixFirst;

                // Run the optimisation
                bool success = optim::lbfgs(x, optFunc, &opt_data, settings);

                // Run the final values
                setAngles(chains, x);
                double finalObj = getObjective(chains);
                if (finalObj < bestObjective) {
                    bestObjective = finalObj;
                }
                numPerformed++;
                std::cout << "sqr= " << finalObj << "  bst= " << bestObjective << "  num= " << numPerformed << "  per=" << itersPerFrame << "  inf= " << chains.size() << " " << numFree << " " << numFree*d << "     \n" << std::flush;

                // Start again with random angles
                for (auto& chain : chains) {
                    std::vector<double> angles = chain.getAngles();
                    for (unsigned long int j=0; j<angles.size(); ++j) {
                        angles[j] = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
                    }
                    if (fixFirst) {
                        angles[0] = 0.0;
                    }
                    chain.setAngles(angles);
                }

                // If we've done enough
                if (numPerformed >= maxIters) {
                    modes.erase(modes.begin());
                    break;
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // Optimise to find local minima, then following the gradient of the minima
        } else if (mode == "minima") {

            // If we haven't set the best minima, assume the current point
            if (bestMinimaAngles.size() == 0) {
                bestMinimaAngles = getAngles(chains);
            }

            // Populate the list of angles to test
            if (anglesToTest.size() == 0) {

                // For now just try each angle in each direction
                for (int i=0; i<100; i++) {

                    // Generate a random vector
                    std::vector<double> randVec;
                    for (unsigned long int j=0; j<bestMinimaAngles.size(); ++j) {
                        randVec.push_back(rand() / double(RAND_MAX) * 2.0 * M_PI);
                    }

                    // Normalise the vector to have a mag of delMinima
                    double norm = 0.0;
                    for (unsigned long int j=0; j<randVec.size(); ++j) {
                        norm += randVec[j] * randVec[j];
                    }
                    norm = sqrt(norm);
                    for (unsigned long int j=0; j<randVec.size(); ++j) {
                        randVec[j] *= delMinima / norm;
                    }

                    // Add the test angles
                    std::vector<double> testAngles = bestMinimaAngles;
                    for (unsigned long int j=0; j<randVec.size(); ++j) {
                        testAngles[j] += randVec[j];
                    }
                    anglesToTest.push_back(testAngles);

                }

                // Set the current angles to the first in the list
                angleIndex = 0;
                setAngles(chains, anglesToTest[angleIndex]);

            }

            // Calculate the objective
            double bestObjective = getObjective(chains);
            double maxDiff = 0.0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                maxDiff = std::max(maxDiff, chains[i].getObjective());
            }
            std::cout << "sqr=" << bestObjective << "  max=" << maxDiff << "  del=" << checkDelta << "  ind=" << angleIndex << "  del=" << delMinima << "  bst=" << bestMinimaValue << "     \r" << std::flush;

            // Stop if we have something close
            if (bestObjective < 1) {
                modes.erase(modes.begin());
                continue;
            }

            // Do a few iterations per display frame
            for (int j=0; j<itersPerFrame; j++) {

                // Keep track of each time we update the best objective
                int numChanges = 0;

                // For each non-fixed chain
                for (unsigned long int j=0; j<chains.size(); ++j) {
                    if (chains[j].isFixed()) {
                        continue;
                    }

                    // Check the effect of moving each direction
                    std::vector<double> angles = chains[j].getAngles();
                    for (unsigned long int k=0; k<angles.size(); ++k) {
                        if (k == 0 && fixFirst) {
                            continue;
                        }

                        // Save the old angle
                        double oldAngle = angles[k];

                        // Check plus checkDelta
                        angles[k] += checkDelta;
                        chains[j].setAngles(angles);
                        double objective = getObjective(chains);
                        if (objective >= bestObjective) {
                            angles[k] = oldAngle;
                        } else {
                            bestObjective = objective;
                            numChanges++;
                            continue;
                        }

                        // Check minus checkDelta
                        angles[k] -= checkDelta;
                        chains[j].setAngles(angles);
                        objective = getObjective(chains);
                        if (objective >= bestObjective) {
                            angles[k] = oldAngle;
                        } else {
                            bestObjective = objective;
                            numChanges++;
                            continue;
                        }

                    }
                    chains[j].setAngles(angles);

                }

                // If no changes, decrease checkDelta
                if (numChanges == 0) {
                    checkDelta *= 0.1;

                    // If it's small enough, stop this with opt
                    if (checkDelta < 1e-4) {
                        checkDelta = 3.0;

                        // If we've found a new minima, start again from there
                        if (bestObjective < bestMinimaValue - 0.1) {
                            bestMinimaValue = bestObjective;
                            bestMinimaAngles = getAngles(chains);
                            delMinima = 10.0;
                            anglesToTest.clear();

                        // Otherwise try a different direction
                        } else {
                            angleIndex++;
                            if (angleIndex >= int(anglesToTest.size())) {
                                delMinima *= 10;
                                anglesToTest.clear();
                            } else {
                                setAngles(chains, anglesToTest[angleIndex]);
                            }
                        }

                        break;

                    }
                }

            }

            // Update the visuals
            if (visual) {
                for (long unsigned int i=0; i<chains.size(); ++i) {
                    chains[i].updatePointsAndLines();
                }
            }

        // If told to output all angles in degrees
        } else if (mode == "output" || mode == "outputDeg") {
            
            // Set the output precision
            std::cout << std::setprecision(10);

            // Output the angles of each chain as a new line
            std::vector<double> uniqueAngles;
            int numUnique = 0;
            int numTotal = 0;
            for (unsigned long int i=0; i<chains.size(); ++i) {
                std::vector<double> angles = chains[i].getAngles();
                for (auto& angle : angles) {
                    while (angle < 0.0) {
                        angle += 2.0 * M_PI;
                    }
                    while (angle >= 2.0 * M_PI) {
                        angle -= 2.0 * M_PI;
                    }
                    if (mode == "outputDeg") {
                        angle *= 180.0 / M_PI;
                    }
                }
                std::pair<int, int> vecIndices = chains[i].getVecIndices();
                std::cout << "chain " << i << " (" << vecIndices.first << "," << vecIndices.second << "): ";
                for (auto& angle : angles) {
                    std::cout << angle << " ";
                }
                std::cout << std::endl;
                for (auto& angle : angles) {
                    numTotal++;
                    bool found = false;
                    for (auto& uniqueAngle : uniqueAngles) {
                        if ((mode == "outputDeg" && std::abs(angle - uniqueAngle) < 0.1) || (mode == "output" && std::abs(angle - uniqueAngle) < 0.001)) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        numUnique++;
                        uniqueAngles.push_back(angle);
                    }
                }
            }
            std::sort(uniqueAngles.begin(), uniqueAngles.end());
            std::cout << std::setprecision(5);
            std::cout << "num angles: " << numTotal << std::endl;
            std::cout << "num unique angles: " << numUnique << std::endl;
            std::cout << "unique angles: ";
            for (auto& angle : uniqueAngles) {
                std::cout << angle << " ";
            }

            // Reset the precision
            std::cout << std::setprecision(5);

            // We're done with this mode
            modes.erase(modes.begin());

        // If told to output the true angles of the vectors
        } else if (mode == "outputTrue" || mode == "outputTrueDeg" || mode == "outputVectors") {
            
            // Set the output precision
            std::cout << std::setprecision(10);

            // Get various numbers
            int numVectors = 0;
            for (unsigned long int i=1; i<N.size(); i++) {
                numVectors += N[i];
            }
            int numAnglesPer = N[0];
            std::vector<std::vector<double>> trueAngles(numVectors, std::vector<double>(numAnglesPer, 0.0));

            // For each angle
            for (int i=0; i<N[0]; i++) {

                // For each chain, we have chain_12 = theta_1 + theta_2
                for (unsigned long int j=0; j<chains.size(); j++) {
                    std::pair<int, int> vecIndices = chains[j].getVecIndices();
                    double angle = chains[j].getAngles()[i];
                    while (angle > 2.0*M_PI) {
                        angle -= 2.0*M_PI;
                    }
                    while (angle < 0.0) {
                        angle += 2.0*M_PI;
                    }
                    if (vecIndices.second == 0) {
                        trueAngles[vecIndices.first][i] = angle;
                    } else if (vecIndices.first == 0) {
                        trueAngles[vecIndices.second][i] = -angle;
                    }
                }

            }

            // Output each angle
            if (mode == "outputTrue") {
                for (int i=0; i<numVectors; i++) {
                    std::cout << "vector " << i << ": ";
                    for (int j=0; j<numAnglesPer; j++) {
                        std::cout << trueAngles[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            } else if (mode == "outputTrueDeg") {
                for (int i=0; i<numVectors; i++) {
                    std::cout << "vector " << i << ": ";
                    for (int j=0; j<numAnglesPer; j++) {
                        std::cout << trueAngles[i][j] * 180.0 / M_PI << " ";
                    }
                    std::cout << std::endl;
                }
            } else if (mode == "outputVectors") {
                std::vector<std::vector<std::complex<double>>> trueVectors(numVectors, std::vector<std::complex<double>>(numAnglesPer));
                for (int i=0; i<numVectors; i++) {
                    std::cout << "vector " << i << ": ";
                    for (int j=0; j<numAnglesPer; j++) {
                        double coeff = 1.0 / std::sqrt(d);
                        trueVectors[i][j] = coeff * std::exp(std::complex<double>(0.0, trueAngles[i][j]));
                        std::cout << trueVectors[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                // Calculate all the inner products (gram matrix)
                Eigen::MatrixXd gramMatrix(numVectors, numVectors);
                for (int i=0; i<numVectors; i++) {
                    for (int j=0; j<numVectors; j++) {
                        std::complex<double> inner = 0.0;
                        for (int k=0; k<numAnglesPer; k++) {
                            inner += std::conj(trueVectors[i][k]) * trueVectors[j][k];
                        }
                        gramMatrix(i, j) = std::abs(inner);
                        if (gramMatrix(i, j) < 1e-10) {
                            gramMatrix(i, j) = 0.0;
                        }
                    }
                }
                std::cout << std::setprecision(5);
                std::cout << "Gram matrix:" << std::endl;
                std::cout << gramMatrix << std::endl;

            }

            // Reset the precision
            std::cout << std::setprecision(5);

            // We're done with this mode
            modes.erase(modes.begin());

		// Otherwise just update if the user changed something
		} else if (somethingChanged) {
			for (long unsigned int i=0; i<chains.size(); ++i) {
				chains[i].update();
			}
		}

		// Clear the window, white background
        if (visual) {
            window.clear(sf::Color::White);
            for (auto& chain : chains) {
                chain.draw(window, font);
            }
            for (unsigned long int i=0; i<buttons.size(); ++i) {
                window.draw(buttons[i]);
                window.draw(buttonTexts[i]);
            }
            window.draw(fpsCounter);
            window.draw(iterCounter);
            window.display();
        }

    }

    std::cout << std::endl;

    return 0;
}

