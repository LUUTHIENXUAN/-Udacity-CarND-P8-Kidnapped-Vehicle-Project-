/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "libkdtree/kdtree++/kdtree.hpp"
#include "particle_filter.h"

using namespace std;

/* Initilization Step
 * TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	     x, y, theta and their uncertainties from GPS) and all weights to 1.
	     Add random Gaussian noise to each particle.
	     NOTE: Consult particle_filter.h for more information about this method (and others in this file).
 */

void ParticleFilter::init(double x, double y, double theta, double std[]) {

	///* Set the number of particles
	num_particles = 1000;
	optimize = true;

	///* Use a non deterministic seed
	random_device rd;
	default_random_engine gen(rd());

	///* Standard deviations for x, y, and theta
	double std_x, std_y, std_theta;
	///* GPS measurement uncertainty [x [m], y [m], theta [rad]]
	std_x     = std[0];
	std_y     = std[1];
	std_theta = std[2];

	if (!is_initialized){

		///* Create a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(x, std_x);

		///* Create a normal (Gaussian) distribution for y
		normal_distribution<double> dist_y(y, std_y);

		///* Create a normal (Gaussian) distribution for theta
		normal_distribution<double> dist_theta(theta, std_theta);

		///* set size of vector particles:
		particles.resize(num_particles);
		///* set size of weights:
		weights.resize(num_particles);

		for ( int i = 0; i < num_particles; ++i) {

			particles[i].id     = i;
			particles[i].x      = dist_x(gen);
			particles[i].y      = dist_y(gen);
			particles[i].theta  = dist_theta(gen);
			particles[i].weight = 1.0;

			weights[i]           = 1;

	    }

		// done initializing
		is_initialized = true;
		return;
    }

}

/* Prediction Step
 * TODO: Add measurements to each particle and add random Gaussian noise.
	     NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	     http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	     http://www.cplusplus.com/reference/random/default_random_engine/
 */

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {

	///* Use a non deterministic seed
	random_device rd;
	mt19937 gen(rd());

	for (auto & p : particles) {

		if(fabs(yaw_rate) < 0.0001){
			///* Updating x, y and the yaw angle when the yaw rate is zero

			p.x     = p.x + velocity * delta_t * cos(p.theta);
			p.y     = p.y + velocity * delta_t * sin(p.theta);
			p.theta = p.theta;
		} else {
			///* Updating x, y and the yaw angle when the yaw rate is not equal to zero

			p.x     = p.x + velocity/yaw_rate * (sin(p.theta + (yaw_rate * delta_t)) - sin(p.theta));
			p.y     = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + (yaw_rate * delta_t)));
			p.theta = p.theta + yaw_rate * delta_t;
		}

		///* Add measurements to each particle and add random Gaussian noise.

		std::normal_distribution<double> noise_x(0, std_pos[0]);
	  std::normal_distribution<double> noise_y(0, std_pos[1]);
	  std::normal_distribution<double> noise_psi(0, std_pos[2]);

		p.x      = p.x + noise_x(gen);
		p.y      = p.y + noise_y(gen);
		p.theta  = p.theta + noise_psi(gen);

	}

}

/* Update Step
 * dataAssociation
 * TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	     observed measurement to this particular landmark.
	     NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
         implement this method and use it as a helper during the updateWeights phase.
 * 1. Transform each observation to the global coordinate system
 * 2. For each transformed observation, find a landmark from the map which is closest to it.
 * 3. Compute the weight for the transformed observation using its closest landmark.
 * 4. Multiply weights of all observations together to get the final weight of the particle.
 */

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
                                     std::vector<LandmarkObs>& observations) {


}

/* updateWights
 * TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	     more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	     NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	     according to the MAP'S coordinate system. You will need to transform between the two systems.
	     Keep in mind that this transformation requires both rotation AND translation (but no scaling).

	     The following is a good resource for the theory:
	     https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	     and the following is a good resource for the actual equation to implement (look at equation 3.33
	     http://planning.cs.uiuc.edu/node99.html
*/

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		                               const std::vector<LandmarkObs> &observations,
		                               const Map &map_landmarks) {

	// clear weights vector
  weights.clear();
	cout << "UpdateWeights Debug" << endl;


	for (auto & p : particles){

		// Restart this particle's weight
		p.weight = 1;

		///* Gather all landmarks inside the sensor_range at one particle
		vector<LandmarkObs> landmarks_inrange;
		landmarks_inrange.reserve(map_landmarks.landmark_list.size());

		for (auto & landmarks : map_landmarks.landmark_list){

			double distance = dist(p.x, p.y, landmarks.x_f, landmarks.y_f);

			if (distance < sensor_range){

				LandmarkObs this_landmark;
				this_landmark.id = landmarks.id_i;
				this_landmark.x  = landmarks.x_f ;
				this_landmark.y  = landmarks.y_f ;

				landmarks_inrange.push_back(this_landmark);

			}
		}


    ///* Transform the car's measurements to the map's coordinate system for one particle
    vector<LandmarkObs> transformed_observations;
    transformed_observations.reserve(observations.size());

		for (auto& obs: observations) {

			///* Transform each observation to the global coordinate system
			LandmarkObs transformed_obs;
			transformed_obs.x  = p.x + cos(p.theta)* obs.x  -  sin(p.theta)* obs.y;
			transformed_obs.y  = p.y + sin(p.theta)* obs.x  +  cos(p.theta)* obs.y;
			transformed_obs.id = obs.id;

			transformed_observations.push_back(transformed_obs);

		}

		///* For each transformed observation, find a closet landmark from the map.
		for (auto& obs: transformed_observations) {

			vector<double> closet_distances;
			closet_distances.reserve(landmarks_inrange.size());

			LandmarkObs predicted_obs;

			for (auto & landmark : landmarks_inrange){

				double this_distance = dist(obs.x, obs.y, landmark.x, landmark.y);

				closet_distances.push_back(this_distance);
				double smallest = *min_element(begin(closet_distances), end(closet_distances));

				if (this_distance == smallest){

					predicted_obs.x  = landmark.x;
					predicted_obs.y  = landmark.y;
					predicted_obs.id = obs.id;

				}
			}

			///* Update weight for this observation
			///* calculate normalization term
			double gauss_norm = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));

			///* calculate exponent
			double exponent   =	((obs.x - predicted_obs.x)*(obs.x - predicted_obs.x))/(2.0 * std_landmark[0] * std_landmark[0]) +
					                ((obs.y - predicted_obs.y)*(obs.y - predicted_obs.y))/(2.0 * std_landmark[1] * std_landmark[1]);

		  ///* calculate weight using normalization terms and exponent
		  double weight     = gauss_norm * exp(-exponent);

		  //* Particle's Final Weight Update
		  p.weight         *= weight;

		}

		weights.push_back(p.weight);

    }


}

void ParticleFilter::updateWeights_optimized(double sensor_range, double std_landmark[],
		                                         const std::vector<LandmarkObs_kd> &observations,
		                                         const Map &map_landmarks) {

	// clear weights vector
  weights.clear();

	for (auto & p : particles){

		KDTree::KDTree<2, LandmarkObs_kd, tac> mapSearch;
    KDTree::KDTree<2, LandmarkObs_kd, tac> map_inrange;
		// Restart this particle's weight
		p.weight = 1;

		///* Gather all landmarks inside the sensor_range at one particle
		for (auto & landmarks : map_landmarks.landmark_list){
			map_inrange.insert(LandmarkObs_kd(landmarks.id_i, landmarks.x_f, landmarks.y_f));
		}

		vector<LandmarkObs_kd> landmarks_inrange;
		LandmarkObs_kd ref;
		ref.id = 0; ref.dim[0]  = p.x; ref.dim[1]  = p.y;
    map_inrange.find_within_range(ref,sensor_range,back_insert_iterator<vector<LandmarkObs_kd> >(landmarks_inrange));

		for (auto & landmark : landmarks_inrange){
			mapSearch.insert(landmark);
		}

		// once you add all your landmarks, you optimise the current tree to search faster
		mapSearch.optimise();

    ///* Transform the car's measurements to the map's coordinate system for one particle
    vector<LandmarkObs_kd> transformed_observations;
    transformed_observations.reserve(observations.size());

		for (auto& obs: observations) {

			///* Transform each observation to the global coordinate system
			LandmarkObs_kd transformed_obs;
			transformed_obs.dim[0]  = p.x + cos(p.theta)* obs.dim[0]  -  sin(p.theta)* obs.dim[1];
			transformed_obs.dim[1]  = p.y + sin(p.theta)* obs.dim[0]  +  cos(p.theta)* obs.dim[1];
			transformed_obs.id = obs.id;

			transformed_observations.push_back(transformed_obs);

		}

		///* For each transformed observation, find a closet landmark from the map.
		for (auto& obs: transformed_observations) {

			// find_nearest will return an iterator that contains a tuple,
			// the first element contains the element that you are looking for.

      LandmarkObs_kd nearest = *mapSearch.find_nearest(obs).first;

			///* Update weight for this observation
			///* calculate normalization term
			double gauss_norm = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));

			///* calculate exponent
			double exponent   =
			((obs.dim[0] - nearest.dim[0])*(obs.dim[0] - nearest.dim[0]))/(2.0 * std_landmark[0] * std_landmark[0]) +
			((obs.dim[1] - nearest.dim[1])*(obs.dim[1] - nearest.dim[1]))/(2.0 * std_landmark[1] * std_landmark[1]);

		  ///* calculate weight using normalization terms and exponent
		  double weight     = gauss_norm * exp(-exponent);

		  //* Particle's Final Weight Update
		  p.weight         *= weight;

		}

		weights.push_back(p.weight);

    }

}

/*
 *TODO: Resample particles with replacement with probability proportional to their weight.
        NOTE: You may find std::discrete_distribution helpful here.
        http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
*/
void ParticleFilter::resample() {

	random_device rd;
    mt19937 gen(rd());

    ///* sample particles based on their weight
    discrete_distribution<> distribution(weights.begin(), weights.end());


	vector<Particle> resample_particles;
	resample_particles.reserve(particles.size());

	for (int i = 0; i < num_particles; i++) {

        Particle sample = particles[distribution(gen)];
        sample.id = i;
        resample_particles.push_back(sample);

  }


    //particles = resample_particles;
		particles = std::move(resample_particles);


}

Particle ParticleFilter::SetAssociations(Particle& particle,
                                         const std::vector<int>& associations,
                                         const std::vector<double>& sense_x,
                                         const std::vector<double>& sense_y){

    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
