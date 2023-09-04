#include <SFML/Graphics.hpp>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>

class Chain {
private:

	// Things needed regarding of visualisation
	int numVectors = 0;
    int dragging = -1;
	double radius = 0;
	int vectorLength;
	double objective = 0.0;
	bool fixFirst = true;
	std::vector<double> angles;
	std::vector<std::pair<double,double>> positions;
	std::pair<double,double> startPosition;
	std::pair<int,int> vecIndices;
	std::vector<std::pair<double, Chain*>> relation;

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

	// Update the points and lines based on the angles
	void updatePointsAndLines() {

		// Update the points based on the angles
		std::pair<double,double> prevPos = startPosition;
		for (int i = 0; i < points.size(); ++i) {
			positions[i] = {prevPos.first + vectorLength * sin(angles[i]), prevPos.second - vectorLength * cos(angles[i])};
			points[i].setPosition(positions[i].first, positions[i].second);
			prevPos = positions[i];
		}

		// Update the lines based on the points
		prevPos = startPosition;
		for (int i = 0; i < lines.size(); ++i) {
			double angle = atan2(positions[i].second - prevPos.second, positions[i].first - prevPos.first);
			lines[i].setPosition(prevPos.first, prevPos.second);
			lines[i].setRotation(angle * 180.0 / M_PI);
			prevPos = positions[i];
		}

	}

	// Set how this chain is defined related to the others
	void setRelation(std::vector<std::pair<double, Chain*>> relation_) {
		relation = relation_;

		// Also change the colour to make it more obvious
		for (int i = 0; i < points.size(); ++i) {
			points[i].setFillColor(sf::Color::Red);
		}

	}

	// Snap the end point to the radius
	void snap() {

		// How close the last point is versus the circle radius
		double finalXMinus2 = 0;
		double finalYMinus2 = 0;
		for (int i = 0; i < angles.size()-2; ++i) {
			finalXMinus2 += vectorLength * sin(angles[i]);
			finalYMinus2 -= vectorLength * cos(angles[i]);
		}
		double finalXAll = finalXMinus2 + vectorLength * sin(angles[angles.size()-2]) + vectorLength * sin(angles[angles.size()-1]);
		double finalYAll = finalYMinus2 - vectorLength * cos(angles[angles.size()-2]) - vectorLength * cos(angles[angles.size()-1]);
		double finalDistMinus2 = sqrt(finalXMinus2 * finalXMinus2 + finalYMinus2 * finalYMinus2)- radius;
		double finalDist = sqrt(finalXAll * finalXAll + finalYAll * finalYAll) - radius;
		double finalAngle = atan2(finalYAll, finalXAll);

		// Move this amount
		move({-finalDist * cos(finalAngle), -finalDist * sin(finalAngle)});
		return;

		// If the final distance is less than the radius, snap it
		if (std::abs(finalDist) < vectorLength / 2.0) {

			// Get the angle between the last two points
			double angle1OG = angles[angles.size()-1];
			double angle2OG = angles[angles.size()-2];
			double angleDiff = angle1OG - angle2OG;
			while (angleDiff < -M_PI) {
				angleDiff += 2*M_PI;
			}
			while (angleDiff > M_PI) {
				angleDiff -= 2*M_PI;
			}
			
			// Do a binary search to find the correct angle
			double lowerBound = -M_PI / 4.0;
			double upperBound = M_PI / 4.0;
			for (int j=0; j<30; j++) {

				// Try the mid point
				double testPoint = (lowerBound + upperBound) / 2.0;

				// Bend in the correct way
				if (angleDiff > 0) {
					angles[angles.size()-1] = angle1OG - testPoint;
					angles[angles.size()-2] = angle2OG + testPoint;
				} else {
					angles[angles.size()-1] = angle1OG + testPoint;
					angles[angles.size()-2] = angle2OG - testPoint;
				}

				// Get the new distance
				double newX = finalXMinus2;
				double newY = finalYMinus2;
				for (int i = angles.size()-2; i < angles.size(); ++i) {
					newX += vectorLength * sin(angles[i]);
					newY -= vectorLength * cos(angles[i]);
				}
				double newDist = sqrt(newX * newX + newY * newY) - radius;

				// Set the bounds
				if (newDist > 0) {
					upperBound = testPoint;
				} else {
					lowerBound = testPoint;
				}

			}

			// Update everything
			update();

		}

	}

	// Getter
	std::pair<int, int> getVecIndices() {
		return vecIndices;
	}

	// Check if the user is dragging any point
	bool isDragging() {
		return dragging != -1;
	}

	// Check if the user can modify this chain
	bool isFixed() {
		return relation.size() != 0;
	}

	// Update the objective
	void updateObjective() {

		// How close the last point is versus the circle radius
		double finalX = 0;
		double finalY = 0;
		for (int i = 0; i < angles.size(); ++i) {
			finalX += vectorLength * sin(angles[i]);
			finalY -= vectorLength * cos(angles[i]);
		}
		double finalDist = sqrt(finalX * finalX + finalY * finalY);
		objective = abs(finalDist - radius);

		// Set the label text
		text.setString(std::to_string(objective));

	}

	// Main constructor
    Chain(int numVectors_, std::pair<double,double> startPosition_, int vectorLength_, int circleRadius_, std::pair<int,int> vecIndices_, bool fixFirst_) {

		// Set the size of the points
		numVectors = numVectors_;
		vectorLength = vectorLength_;
		startPosition = startPosition_;
		vecIndices = vecIndices_;
		radius = circleRadius_;
		fixFirst = fixFirst_;

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
			sf::RectangleShape line(sf::Vector2f(vectorLength, lineThickness));
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

	// Given a delta vector, move a given point
	void move(std::pair<double,double> delta, int toDrag=-1) {

		// If told to drag the last point
		if (toDrag == -1) {
			toDrag = points.size() - 1;
		}

		// The points after the dragged point
		for (int i=toDrag; i<points.size(); ++i) {

			// Simple move with the point
			positions[i] = std::pair<double,double>(positions[i].first + delta.first, positions[i].second + delta.second);

		}

		// The points before the dragged point
		for (int i=toDrag-1; i>=0; --i) {

			// Get the angle between this point at the one after it
			double angles = atan2(positions[i + 1].second - positions[i].second, positions[i + 1].first - positions[i].first);
			
			// Move with that angle to the right distance
			positions[i] = std::pair<double,double>(positions[i + 1].first - vectorLength * cos(angles), positions[i + 1].second - vectorLength * sin(angles));

		}

		// Set all the angles
		std::pair<double,double> prevPos = startPosition;
		for (int i=0; i<points.size(); ++i) {
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
			for (int i = 0; i < points.size(); i++) {
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
	void update() {

		// If it's fixed, update the angles
		if (relation.size() != 0) {

			// For each element in the chain
			for (int i = 0; i < points.size(); ++i) {

				// Get the sum of the angles from the other chains
				double angle = 0;
				for (int j = 0; j < relation.size(); ++j) {
					angle += relation[j].first * relation[j].second->angles[i];
				}

				// Set the angle to this point
				angles[i] = angle;

			}

		}

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
		for (int i = 0; i < points.size(); ++i) {
            window.draw(points[i]);
        }
		text.setFont(font);
		window.draw(text);
    }

};

// Entry point
int main(int argc, char* argv[]) {

	// Settings
	int d = 3;
	std::vector<int> N = {d, 1, 1, 1};
	bool fixFirst = true;
	bool anneal = false;
	double startTemp = 1.0;
	int steps = 1000;

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

			// Require that the first basis size is the same as the dimension
			d = N[0];

		// If asked to display the help
		} else if (std::string(argv[i]) == "-h") {
			std::cout << "Usage: " << argv[0] << " [-N <basis sizes>] [-a <temp> <steps>] [-1]" << std::endl;
			return 0;

		// If told not to fix the first angle
		} else if (std::string(argv[i]) == "-1") {
			fixFirst = false;

		// If told to anneal
		} else if (std::string(argv[i]) == "-a") {
			anneal = true;
			startTemp = std::stof(argv[i + 1]);
			steps = std::stoi(argv[i + 2]);

		}

	}

	// Create the main window
	sf::ContextSettings settings;
    settings.antialiasingLevel = 3.0;
	int windowWidth = 1280;
	int windowHeight = 720;
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Visual MUBs", sf::Style::Default, settings);
	window.setFramerateLimit(60);

	// Load the font
	sf::Font font;
	if (!font.loadFromFile("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")) {
		std::cout << "Error loading font" << std::endl;
	}

	// Create each of the chains representing the MUBs
	double scaling = 100.0;
	double spacing = scaling * 3.0;
	int gridX = 0;
	int gridY = 0;
	double minX = 0;
	double minY = 0;
	double maxX = 0;
	double maxY = 0;
	int NSoFarI = 0;
	int NSoFarJ = 0;
	std::vector<Chain> chains;
	for (int i = 1; i < N.size(); ++i) {
		int NSoFarI = 0;
		for (int m = 1; m < i; ++m) {
			NSoFarI += N[m];
		}
		for (int j = i; j < N.size(); ++j) {
			int NSoFarJ = 0;
			for (int m = 1; m < j; ++m) {
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
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, 0.0, {gridX, gridY}, fixFirst));

					// Mutually unbiasedness
					} else {
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, scaling, {gridX, gridY}, fixFirst));
					}

				}
			}
		}
	}

	// Create an empty eigen matrix
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(std::pow(chains.size(), 2), chains.size());
	int nextInd = 0;

	// Add the linear constraints
	for (int i = 0; i < chains.size(); ++i) {
		std::pair<int,int> vecIndices1 = chains[i].getVecIndices();	
		for (int j = i+1; j < chains.size(); ++j) {
			std::pair<int,int> vecIndices2 = chains[j].getVecIndices();

			// theta_abi + theta_cai = theta_cbi
			if (vecIndices1.first == vecIndices2.second) {

				// Find the other
				std::pair<int,int> lookingFor = std::make_pair(vecIndices2.first, vecIndices1.second);
				int otherInd = -1;
				for (int k = 0; k < chains.size(); ++k) {
					if (chains[k].getVecIndices() == lookingFor) {
						otherInd = k;
						break;
					}
				}

				// Add to the matrix
				if (otherInd >= 0) {
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

	// Flip the matrix horizontally
	A = A.rowwise().reverse().eval();

	// Row reduce the matrix A using Eigen
	Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
	Eigen::MatrixXd reducedA = lu.matrixLU().triangularView<Eigen::Upper>();
	int rank = lu.rank();
	reducedA.conservativeResize(rank, Eigen::NoChange);

	// Make sure each of the diagonals is 1
	for (int i = 0; i < reducedA.rows(); ++i) {
		reducedA.row(i) /= reducedA(i,i);
	}

	// Make sure each of the diagonals is the only one in that column
	for (int i=0; i<reducedA.rows(); ++i) {
		for (int j=i+1; j<reducedA.rows(); ++j) {
			if (abs(reducedA(i,j)) > 1e-10) {
				reducedA.row(i) -= reducedA(i,j)*reducedA.row(j);
			}
		}
	}

	// For each row of the reduced matrix, get the first one
	for (int i = 0; i < reducedA.rows(); ++i) {
		for (int j = 0; j < reducedA.cols(); ++j) {
			if (reducedA(i,j) != 0) {
				std::vector<std::pair<double, Chain*>> terms;
				std::cout << j << " ";
				for (int k = j+1; k < reducedA.cols(); ++k) {
					if (std::abs(reducedA(i,k)) > 1e-10) {
						terms.push_back(std::make_pair(reducedA(i,k), &chains[k]));
						std::cout << k << " ";
					}
				}
				std::cout << std::endl;
				chains[j].setRelation(terms);
				break;
			}
		}
	}

	// Vars used for event management
	bool draggingBackground = false;
	sf::Vector2i lastMousePos;
	sf::View currentView = window.getView();
	sf::Vector2f lastViewPos;
	double zoomLevel = 1.0;
	double currentTemp = startTemp;
	double deltaTemp = startTemp / steps;
	double prevObjective = 100000000.0;

	// Set the initial view so we can see everything
	currentView.setCenter((minX+maxX+windowWidth+scaling)/2.0, (minY+maxY+windowHeight+scaling)/2.0);
	window.setView(currentView);

	// Start the main loop
    while (window.isOpen()) {

		// Process events
        sf::Event event;
		bool somethingChanged = false;
        while (window.pollEvent(event)) {

			// Each chain should process events first
            for (auto& chain : chains) {
				bool val = chain.handleEvent(event, window);
				somethingChanged = somethingChanged || val;
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

		// If told to anneal TODO
		if (anneal && currentTemp > 0.0) {

			// Pick a random chain that is not fixed
			int chainIndex = rand() % chains.size();
			while (chains[chainIndex].isFixed()) {
				chainIndex = rand() % chains.size();
			}

			// Pick a random distance proportional to the temperature
			double distance = currentTemp * (rand() / (double)RAND_MAX);
			double angle = 2.0 * M_PI * (rand() / (double)RAND_MAX);

			// Move the chain that amount
			std::pair<double, double> movement = {distance * cos(angle), distance * sin(angle)};
			chains[chainIndex].move(movement);

			// See if the overall objective function is better
			double overallObjective = 0.0;
			for (auto& chain : chains) {
				chain.update();
				overallObjective += std::pow(chain.getObjective(), 2);
			}

			// Accept with Boltzmann probability
			if (overallObjective < prevObjective || (rand() / (double)RAND_MAX) < std::exp((prevObjective - overallObjective) / currentTemp)) {
				prevObjective = overallObjective;
			} else {
				chains[chainIndex].move({-movement.first, -movement.second});
			}

			// Lower the temperature
			currentTemp -= deltaTemp;

			std::cout << overallObjective << " " << currentTemp << std::endl;

		// Apply the linear constraints to each chain
		} else if (somethingChanged) {
			for (int i = 0; i < chains.size(); ++i) {
				chains[i].update();
			}
		}

		// Clear the window, white background
        window.clear(sf::Color::White);
		for (auto& chain : chains) {
			chain.draw(window, font);
		}
        window.display();

    }

    return 0;
}

