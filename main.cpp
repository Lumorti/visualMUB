#include <SFML/Graphics.hpp>
#include <math.h>
#include <iostream>

class Chain {
private:
    std::vector<sf::CircleShape> points;
    std::vector<sf::RectangleShape> lines;
    int dragging = -1;
	int pointSize = 5;
	int lineThickness = 2;
	sf::Vector2f startPosition;
	sf::CircleShape radiusCircle;
	int vectorLength;
	sf::Text label;

public:

	// Update the lines based on the points
	void updateLines() {
		for (int i = 0; i < lines.size(); ++i) {
			float angle = atan2(points[i + 1].getPosition().y - points[i].getPosition().y, points[i + 1].getPosition().x - points[i].getPosition().x);
			lines[i].setPosition(points[i].getPosition().x + pointSize, points[i].getPosition().y + pointSize);
			lines[i].setRotation(angle * 180 / 3.14159265358979323846);
		}
	}

	// Check if the user is dragging any point
	bool isDragging() {
		return dragging != -1;
	}

	// Main constructor
    Chain(int numVectors, sf::Vector2f startPosition_, int vectorLength_, int circleRadius_, std::string name_) {

		// Set the size of the points
		numVectors = numVectors+1;
		vectorLength = vectorLength_;
		startPosition = startPosition_;

		// Create the points
        for (int i = 0; i < numVectors; ++i) {
            sf::CircleShape point(pointSize);
            point.setFillColor(sf::Color::Blue);
            point.setPosition(startPosition.x + i * vectorLength, startPosition.y);
            points.push_back(point);
        }

		// Create the lines
		for (int i = 0; i < numVectors - 1; ++i) {
			sf::RectangleShape line(sf::Vector2f(vectorLength, lineThickness));
			line.setFillColor(sf::Color::Black);
			line.setPosition(startPosition.x + i * vectorLength, startPosition.y);
			lines.push_back(line);
		}

		// Create the radius circle
		radiusCircle.setRadius(circleRadius_);
		radiusCircle.setOrigin(circleRadius_, circleRadius_);
		radiusCircle.setPosition(startPosition.x+pointSize, startPosition.y+pointSize);
		radiusCircle.setFillColor(sf::Color::Transparent);
		radiusCircle.setOutlineColor(sf::Color::Red);
		radiusCircle.setOutlineThickness(lineThickness);

		// Create the title
		label.setString(name_);
		label.setCharacterSize(20);
		label.setFillColor(sf::Color::Black);
		label.setPosition(startPosition.x - name_.size()*5, startPosition.y - circleRadius_ - 30);

		// Set all the line positions
		updateLines();

    }

	// When an event happens, the chain should handle it
    void handleEvent(const sf::Event& event, sf::RenderWindow& window) {

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

		// When the mouse button is released, the point being dragged should stop
		} else if (event.type == sf::Event::MouseButtonReleased) {
			if (event.mouseButton.button == sf::Mouse::Left) {
				dragging = -1;
			}

		// When the mouse is moved, the point being dragged should follow it
		} else if (event.type == sf::Event::MouseMoved) {

			// If there is a point being dragged
			if (dragging != -1) {

				// Convert to coords
				sf::Vector2f mousePosition = window.mapPixelToCoords(sf::Vector2i(event.mouseMove.x, event.mouseMove.y));

				// Move that point to the mouse, talking into account the size
				double deltaX = mousePosition.x - points[dragging].getPosition().x - pointSize;
				double deltaY = mousePosition.y - points[dragging].getPosition().y - pointSize;
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
					//points[i].setPosition(points[i - 1].getPosition().x - vectorLength * cos(angle), points[i - 1].getPosition().y - vectorLength * sin(angle));
					points[i].move(deltaX, deltaY);

				}

				// Move all points so that point 0 is at the start position
				sf::Vector2f deltaToStart = startPosition - points[0].getPosition();
				for (int i = 0; i < points.size(); ++i) {
					points[i].move(deltaToStart.x, deltaToStart.y);
				}

				// Set all the line positions
				updateLines();

			}

        }

    }

	// Draw the points, lines, circle and labels
    void draw(sf::RenderWindow& window, sf::Font& font) {
		window.draw(radiusCircle);
		label.setFont(font);
		window.draw(label);
		for (const auto& line : lines) {
			window.draw(line);
		}
		for (int i = 1; i < points.size(); ++i) {
            window.draw(points[i]);
        }
    }

};

// Entry point
int main() {

	// Settings
	int d = 3;
	std::vector<int> N = {d, 1, 1, 1};
	int gridWidth = 2;

	// Create the main window
	sf::ContextSettings settings;
    settings.antialiasingLevel = 3.0;
    sf::RenderWindow window(sf::VideoMode(800, 600), "Visual MUBs", sf::Style::Close, settings);
	window.setFramerateLimit(60);

	// Load the font
	sf::Font font;
	if (!font.loadFromFile("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")) {
		std::cout << "Error loading font" << std::endl;
	}

	// Create each of the chains representing the MUBs
	double scaling = 100.0;
	int gridX = 0;
	int gridY = 0;
	std::vector<Chain> chains;
	for (int i = 1; i < N.size(); ++i) {
		for (int j = i; j < N.size(); ++j) {
			for (int k = 0; k < N[i]; ++k) {
				for (int l = 0; l < N[j]; ++l) {

					// We already enforce normality
					if (i == j && k == l) {
						continue;
					}

					// The location of the chain
					float currentX = 200 + gridX*(scaling*2.5);
					float currentY = 200 + gridY*(scaling*2.5);
						
					// Orthogonality
					if (i == j) {
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, 0.0, "ortho"));

					// Mutually unbiasedness
					} else {
						chains.push_back(Chain(d, {currentX, currentY}, sqrt(d)*scaling/d, scaling, "mub"));
					
					}

					// Update the grid coords
					gridX += 1;
					if (gridX >= gridWidth) {
						gridX = 0;
						gridY += 1;
					}

				}
			}
		}
	}

	// Vars used for event management
	bool draggingBackground = false;
	sf::Vector2i lastMousePos;
	sf::View currentView = window.getView();
	sf::Vector2f lastViewPos;
	double zoomLevel = 1.0;

	// Start the main loop
    while (window.isOpen()) {

		// Process events
        sf::Event event;
        while (window.pollEvent(event)) {

			// Each chain should process events first
            for (auto& chain : chains) {
				chain.handleEvent(event, window);
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

		// Clear the window, white background
        window.clear(sf::Color::White);
		for (auto& chain : chains) {
			chain.draw(window, font);
		}
        window.display();

    }

    return 0;
}

