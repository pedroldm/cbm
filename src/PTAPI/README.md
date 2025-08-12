<h1>PT Application Programming Interface</h1>
<h2>Introduction</h2>
<p>The proposed API offers a simple and efficient connection interface, implemented in C++, making it flexible for use across different platforms. The API usage involves three main steps: solution encoding and decoding, coding problem-specific components (evaluation function and movements), and connecting to the PT core implementation.</p>

<h2>Step 1: Solution Encoding</h2>
<p>In this step, a custom C++ structure is created to encode the solution. This structure extends the "solution" structure from the Problem.h file and should contain variables intended for encoding a solution to the problem. Below is an example structure for a permutation problem:</p>

<pre>
struct solPermutation : public solution {
    std::vector<int> sol;
};
</pre>

<h2>Step 2: Problem-Specific Components</h2>
<p>A C++ class is constructed to identify the problem-specific components. This class must implement the Problem interface and include functions for constructing initial solutions, generating neighboring solutions, and evaluating solutions.</p>

<pre>
class Permutation : public Problem {
public:
    Permutation(std::string filename);
    solPermutation construction();
    solPermutation neighbor(solPermutation sol);
    double evaluate(solPermutation sol);
};
</pre>

<h2>Step 3: Connecting to PT Core</h2>
<p>In this step, an object of the previously created class is instantiated, and an object of type PT is created with the desired parameters. Finally, the PT algorithm is started, and the best solution found is returned.</p>

<pre>
#include "../include/PT.h"
#include "../include/Permutation.h"

int main(int argc, char* argv[]) {
    Permutation* prob = new Permutation("F1001.txt");
    PT&lt;solPermutation&gt; algo(0.01, 1.0, 20, 200, 1000, 1, 3, 200);
    solPermutation sol = algo.start(20, prob);
    std::cout << sol.evalSol;
}
</pre>

<h2>Explanation of PT Class Parameters and start Function</h2>

<p>The parameters of the <code>PT</code> class and the <code>start</code> function are as follows:</p>

<ol>
  <li><strong>Initial temperature (float):</strong> Represents the starting temperature value used in the PT process.</li>
  <li><strong>Final temperature (float):</strong> Denotes the final temperature value used in the PT process.</li>
  <li><strong>Number of replicas (int):</strong> Specifies the number of replicas used in the parallel tempering algorithm.</li>
  <li><strong>Length of the homogeneous Markov chain (int):</strong> Indicates the length of the Markov chain used in each replica during the simulation.</li>
  <li><strong>Number of exchanges between temperatures (int):</strong> Defines the number of exchanges attempted between adjacent temperatures in the parallel tempering algorithm.</li>
  <li><strong>Initial temperature distribution (int):</strong> Specifies the method used to distribute initial temperatures across replicas in the parallel tempering algorithm. This parameter can take the following values:
    <ul>
      <li>1 (Linear): Initial temperatures are distributed linearly from the starting temperature to the final temperature.</li>
      <li>2 (Inverse Linear): Initial temperatures are distributed inversely linearly from the final temperature to the starting temperature.</li>
      <li>3 (Exponential): Initial temperatures are distributed exponentially between the starting temperature and the final temperature.</li>
      <li>4 (Geometric Progression): Initial temperatures are distributed in a geometric progression from the starting temperature to the final temperature.</li>
    </ul>
  </li>
  <li><strong>Automatic temperature update (int):</strong> Specifies the method used to automatically update temperatures during the PT process. This parameter can take the following values:
    <ul>
      <li>0 (Disabled): Automatic temperature updating is disabled.</li>
      <li>1 (23%): Temperature is updated to a acceptance rate of 23% .</li>
      <li>2 (Equal Acceptance Rates): Temperature is updated based on equal acceptance rates in adjacent temperatures.</li>
      <li>3 (Feedback Optimized): Temperature is updated using feedback-optimized methods.</li>
    </ul>
  </li>
  <li><strong>Temperature update rate (int):</strong> Represents the rate at which the temperature is updated during the PT process.</li>
</ol>

<p>The <code>start</code> function initializes the execution of the <code>PT</code> algorithm, taking the following parameters:</p>

<ol>
  <li><strong>Number of threads (int):</strong> Specifies the number of threads to be created and used in the parallel tempering algorithm.</li>
  <li><strong>Pointer to the problem object (Problem*):</strong> Points to the object referring to the problem created by the user, which implements the necessary functions for problem-dependent components such as construction, neighbor generation, and evaluation.</li>
</ol>


<h2>Folder Structure</h2>
<pre>
📦Project
 ┣ 📂include
 ┃ ┗ 📜API core files
 ┣ 📂instances
 ┃ ┗ 📜instance files for the examples
 ┣ 📂example
 ┃ ┗ 📜mainSSP.cpp
 ┃ ┗ 📜mainIPMTC.cpp
 ┃ ┗ 📜mainRCPMS.cpp
</pre>

<p>Within the "example" folder, you will find three sample files demonstrating how to use the API for the SSP, IPMTC, and RCPMS problems. Each file contains code examples illustrating how to compile and execute the algorithm to solve each specific problem. These examples are provided to assist users in understanding how to integrate the API into their own projects and utilize its functionalities to solve various optimization problems.</p>

<h2>License</h2>
<p>The API is available for download under the Creative Commons BY-NC license (CC BY-NC) at: <a href="https://creativecommons.org/licenses/">https://creativecommons.org/licenses/</a>. It can be used and adapted for non-commercial purposes, with proper credit given.</p>

