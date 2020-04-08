# COVID-19 analysis based on a SEIR model

This project consists of a disease model for analyzing the spread of COVID-19. 

The project intuition and parameters estimations were retrieved from the article of [Christian Hubbs](https://towardsdatascience.com/social-distancing-to-slow-the-coronavirus-768292f04296). For a more detailed explanation of the modeling process, please read Christian's article.

# Model

The project uses a SEIR model. The SEIR model the flow of the population between 4 states:
- Susceptible
- Exposed
- Infected
- Recovered/Dead

The model assumes that the popualtion has a fixed size and also that its distributed between only these 4 states. Each of these variables represents the number of people in these groups.
The model is governed by 4 differential equations:

![eq1](https://latex.codecogs.com/gif.latex?S%28t%20&plus;%20%5CDelta%20t%29%20%3D%202)
