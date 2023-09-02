#include <SFML/Graphics.hpp>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>

class Chain {
private:
    std::vector<sf::CircleShape> points;
    std::vector<sf::RectangleShape> lines;
	int numVectors = 0;
    int dragging = -1;
	int pointSize = 5;
	int lineThickness = 5;
	double radius = 0;
	sf::Vector2f startPosition;
	sf::CircleShape radiusCircle;
	int vectorLength;
	std::vector<double> angles;
	std::pair<int,int> vecIndices;
	std::vector<std::pair<double, Chain*>> relation;
	double objective = 0.0;
	sf::Text text;

public:

	// Update the points and lines based on the angles
	void updatePointsAndLines() {

		// Update the points based on the angles
		sf::Vector2f prevPos = startPosition;
		for (int i = 0; i < points.size(); ++i) {
			points[i].setPosition(prevPos.x + vectorLength * sin(angles[i]), prevPos.y - vectorLength * cos(angles[i]));
			prevPos = points[i].getPosition();
		}

		// Update the lines based on the points
		prevPos = startPosition;
		for (int i = 0; i < lines.size(); ++i) {
			float angle = atan2(points[i].getPosition().y - prevPos.y, points[i].getPosition().x - prevPos.x);
			lines[i].setPosition(prevPos.x + pointSize, prevPos.y + pointSize);
			lines[i].setRotation(angle * 180.0 / M_PI);
			prevPos = points[i].getPosition();
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

	// Snap the end point to the radius TODO
	void snap() {

		// How close the last point is versus the circle radius
		double finalX = 0;
		double finalY = 0;
		for (int i = 0; i < angles.size()-1; ++i) {
			finalX += vectorLength * sin(angles[i]);
			finalY -= vectorLength * cos(angles[i]);
		}
		double finalDist = sqrt(finalX * finalX + finalY * finalY) - radius;

		// If the final distance is less than the radius, snap it
		if (std::abs(finalDist) < vectorLength) {

			// Spin around until we change sign
			int numChecks = 20;
			double minAngle = 0;
			double maxAngle = 2.0 * M_PI;

			// Keep repeating (each repeat is an order of mag, more or less)
			for (int j=0; j<5; j++) {

				// Search for a sign change in the distance verus the radius
				std::vector<std::pair<double,double>> newAngleRanges;
				double prevAngle = maxAngle;
				double prevDist = sqrt(std::pow(finalX + vectorLength * sin(maxAngle), 2) + std::pow(finalY - vectorLength * cos(maxAngle), 2)) - radius;
				for (int i=0; i<numChecks; i++) {
					double testAngle = minAngle + double(i) * (maxAngle-minAngle) / double(numChecks-1);
					double testDist = sqrt(std::pow(finalX + vectorLength * sin(testAngle), 2) + std::pow(finalY - vectorLength * cos(testAngle), 2)) - radius;
					if (testDist*prevDist < 0 && (i > 0 || j == 0)) {
						newAngleRanges.push_back(std::make_pair(std::min(prevAngle, testAngle), std::max(prevAngle, testAngle)));
					}
					prevDist = testDist;
					prevAngle = testAngle;
				}

				// Check which of the angle ranges is closest to the current angle
				double bestAngleDistance = 1000000;
				int bestRange = -1;
				for (int i=0; i<newAngleRanges.size(); i++) {
					double averageAngle = (newAngleRanges[i].first + newAngleRanges[i].second) / 2.0;
					double angleDiff = angles[angles.size()-1] - averageAngle;
					double angleDistance = std::min(std::abs(angleDiff), 2.0 * M_PI - std::abs(angleDiff));
					if (angleDistance < bestAngleDistance) {
						bestAngleDistance = angleDistance;
						bestRange = i;
					}
				}

				// Update the ranges
				std::cout << "Best range: " << newAngleRanges[bestRange].first << " " << newAngleRanges[bestRange].second << std::endl;
				minAngle = newAngleRanges[bestRange].first;
				maxAngle = newAngleRanges[bestRange].second;

			}

			// Set the angle to the average of the range
			angles[angles.size()-1] = (minAngle + maxAngle) / 2.0;

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
    Chain(int numVectors_, sf::Vector2f startPosition_, int vectorLength_, int circleRadius_, std::pair<int,int> vecIndices_) {

		// Set the size of the points
		numVectors = numVectors_;
		vectorLength = vectorLength_;
		startPosition = startPosition_;
		vecIndices = vecIndices_;
		radius = circleRadius_;

		// For each vector
        for (int i = 0; i < numVectors; ++i) {

			// Create the point
            sf::CircleShape point(pointSize);
            point.setFillColor(sf::Color::Blue);
            points.push_back(point);

			// Create the angle
			angles.push_back(0);

			// Create the line
			sf::RectangleShape line(sf::Vector2f(vectorLength, lineThickness));
			line.setFillColor(sf::Color::Black);
			lines.push_back(line);

        }

		// Create the radius circle
		radiusCircle.setRadius(radius);
		radiusCircle.setOrigin(radius, radius);
		radiusCircle.setPosition(startPosition.x+pointSize, startPosition.y+pointSize);
		radiusCircle.setFillColor(sf::Color::Transparent);
		radiusCircle.setOutlineColor(sf::Color::Red);
		radiusCircle.setOutlineThickness(lineThickness);

		// Create the label
		text.setString(std::to_string(objective));
		text.setCharacterSize(20);
		text.setFillColor(sf::Color::Black);
		text.setPosition(startPosition.x, startPosition.y - radius - 20);

		// Set all the line positions
		update();

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
				float deltaX = mousePosition.x - points[dragging].getPosition().x - pointSize;
				float deltaY = mousePosition.y - points[dragging].getPosition().y - pointSize;
				points[dragging].move(deltaX, deltaY);

				// The points before the dragged point
				for (int i = dragging - 1; i >= 0; --i) {

					// Get the angle between this point at the one after it
					float angle = atan2(points[i + 1].getPosition().y - points[i].getPosition().y, points[i + 1].getPosition().x - points[i].getPosition().x);
					
					// Move with that angle to the right distance
					points[i].setPosition(points[i + 1].getPosition().x - vectorLength * cos(angle), points[i + 1].getPosition().y - vectorLength * sin(angle));

				}

				// The points after the dragged point
				for (int i = dragging + 1; i < points.size(); ++i) {

					// Get the angle between this point at the one before it
					float angle = atan2(points[i - 1].getPosition().y - points[i].getPosition().y, points[i - 1].getPosition().x - points[i].getPosition().x);

					// Move with that angle to the right distance
					points[i].move(deltaX, deltaY);

				}

				// Set all the angles
				sf::Vector2f prevPos = startPosition;
				for (int i = 0; i < points.size(); ++i) {
					angles[i] = M_PI/2.0+atan2(points[i].getPosition().y-prevPos.y, points[i].getPosition().x-prevPos.x);
					prevPos = points[i].getPosition();
				}

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
				float angle = 0;
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
			std::cout << "Usage: " << argv[0] << " [-N <basis sizes>]" << std::endl;
			return 0;

		}

	}

	// Create the main window
	sf::ContextSettings settings;
    settings.antialiasingLevel = 3.0;
	int windowWidth = 1280;
	int windowHeight = 720;
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Visual MUBs", sf::Style::Close, settings);
	window.setFramerateLimit(60);

	// Load the font
	sf::Font font;
	if (!font.loadFromFile("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")) {
		std::cout << "Error loading font" << std::endl;
	}

	// Create each of the chains representing the MUBs
	float scaling = 100.0;
	float spacing = scaling * 3.0;
	int gridX = 0;
	int gridY = 0;
	float minX = 0;
	float minY = 0;
	float maxX = 0;
	float maxY = 0;
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
					float currentX = gridX*spacing;
					float currentY = gridY*spacing;

					// Only the upper triangle
					if (gridX <= gridY) {
						continue;
					}

					// Update the min and max values
					minX = std::min(minX, currentX);
					minY = std::min(minY, currentY);
					maxX = std::max(maxX, currentX);
					maxY = std::max(maxY, currentY);
						
					// Orthogonality TODO putting first causes issues
					if (i == j) {
						//chains.insert(chains.begin(), Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, 0.0, {gridX, gridY}));
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, 0.0, {gridX, gridY}));

					// Mutually unbiasedness
					} else {
						//chains.insert(chains.begin(), Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, scaling, {gridX, gridY}));
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, scaling, {gridX, gridY}));
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

	std::cout << A << std::endl;

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
	float zoomLevel = 1.0;

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

		// Apply the linear constraints to each chain
		if (somethingChanged) {
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

