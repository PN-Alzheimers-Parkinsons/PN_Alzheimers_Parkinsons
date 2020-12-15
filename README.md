# Brandon Maddy Github Repository

## Explanation of code
Note: This was our first time writing code in Python and there are therefore aspects which will not have been the most efficient way in which we could have written it. We hope this page will explain our logic and how the different aspects of the code work.

### Implementation of stochasticity
Our model creates stochasticity in the number of tokens transferred at each transition. This is achieved by drawing a number from a distribution, with a mean equal to the rate function and a standard deviation which can be altered. The distribution_type attribute of the add_transition function is a list which takes a number of parameters relating to the random distribution. The first element defines the distribution with with the number is drawn from; in this case we have used a Guassian distribution. The second element is the standard deviation. The third element is the rate function. This function has been determined from the literature, and takes the number of tokens in the relevent place at  to evaluate token transfer for each time step. The fourth element in the list is the firing condition. This condition must return True for a transition to 'fire' and token transfer to occur. 

We have also implemented a calcium switch into our model. This was done by creating artificial places which transfer a set number of tokens at each time step until a threshold is reached, at which point all tokens are then transferred immediately to the next place.  This creates a periodicity in calcium levels, which simulates cycles of calcium concentration in cells. 

## Link to Shared Google Doc of Road Map
https://docs.google.com/document/d/1_yZp0PqxzRHkoOPmO9z-ACW-awMI4X-gydiXJsOXBC4/edit
