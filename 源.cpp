#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536

using namespace std;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	double fitness;
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	double pBestFitness;

	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
		this->pBestFitness = best_fitness;
	}
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c1, double c2, double wmax, double wmin, double dt, double percent)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c1 = c1;
		this->c2 = c2;
		this->wmax = wmax;
		this->wmin = wmin;
		this->dt = dt;
		this->percent = percent;
		particles.resize(m);
		gBest.resize(dim);
	}


	double fitnessFunction(Particle particle)
	{
		double result = 0.0;
		for (int i = 0; i < dim; i++)
		{
			result += pow(particle.position[i], 2);
		}

		/*
		double temp1 = 0.0;
		double temp2 = 0.0;
		for (int i = 0; i < dim; i++)
		{
			double x = particle.position[i];
			temp1 += pow(x, 2);
			temp2 += cos(2 * pi*x);
		}

		result = -20 * exp(-0.2*sqrt(temp1 / dim)) + 20 + E - exp(temp2 / dim);
		*/
		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].position[j] = particles[i].position[j];

		}
		particles[i].fitness = fitnessFunction(particles[i]);
		particles[i].pBestFitness = fitnessFunction(particles[i]);
	}

	void initialAllParticles()
	{
		initialParticles(0);
		gBest_fitness = particles[0].pBestFitness;
		for (int i = 0; i < dim; i++)
			gBest[i] = particles[0].pBest[i];

		for (int i = 1; i < m; i++)
		{
			initialParticles(i);
			if (particles[i].pBestFitness < gBest_fitness)
			{
				gBest_fitness = particles[i].pBestFitness;
				for (int j = 0; j < dim; j++)
				{
					gBest[j] = particles[i].pBest[j];
				}
			}
		}
	}

	void inertiaWeight()
	{
		//w = randDouble(0.4, 0.6);
		//double t = T / ((double)Tmax);
		//w = wmin + (wmax - wmin)*(t - 1)*(t - 1);
		w = 0.9;
	}

	void updateParticle(int i)
	{
		for (int j = 0; j < dim; j++)
		{
			double last_position = particles[i].position[j];
			double range = percent * (max - min);

			particles[i].velocity[j] = w * particles[i].velocity[j] +
				c1 * randDouble(0, 1) * (particles[i].pBest[j] - particles[i].position[j])
				+ c2 * randDouble(0, 1) * (gBest[j] - particles[i].position[j]);
			particles[i].position[j] += dt * particles[i].velocity[j];

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = randDouble(-range, range);

			if (particles[i].velocity[j] < -range)
			{
				particles[i].velocity[j] = randDouble(-range,range);
			}

			if (particles[i].position[j] > max)
			{
				double thre = randDouble(0,1);
				if (last_position == max)
				{
					particles[i].position[j] = randDouble(min, max);
				}
				else if (thre < 0.5)
				{
					particles[i].position[j] = max - (max - last_position) * randDouble(0,1);
				}
				else
				{
				particles[i].position[j] = max;
				}
			}
			if (particles[i].position[j] < min)
			{
				double thre = randDouble(0,1);
				if (last_position == min)
				{
					particles[i].position[j] = randDouble(min, max);
				}
				else if (thre < 0.5)
				{
					particles[i].position[j] = min + (last_position - min) * randDouble(0, 1);
				}
				else
				
				particles[i].position[j] = min;

			}



		}
		particles[i].fitness = fitnessFunction(particles[i]);
		if (particles[i].fitness < particles[i].pBestFitness)
		{
			particles[i].pBestFitness = particles[i].fitness;
			for (int j = 0; j < dim; j++)
			{
				particles[i].pBest[j] = particles[i].position[j];
			}
		}

	}


	void updateAllParticles()
	{
		inertiaWeight();
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
			if (particles[i].pBestFitness < gBest_fitness)
			{
				gBest_fitness = particles[i].pBestFitness;
				for (int j = 0; j < dim; j++)
				{
					gBest[j] = particles[i].pBest[j];
				}
			}
		}
		T++;
	}

	double getFitness()
	{
		return gBest_fitness;
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c1;
	double c2;
	double wmax;
	double wmin;

	double dt;//时间步长
	double percent;

	vector<double>gBest;
	double gBest_fitness;

	vector<Particle> particles;


};

int main()
{

	int dim = 30;
	int m = 20;
	int Tmax = 2000;
	double max = 100;
	double min = -100;
	double c1 = 2;
	double c2 = 2;
	double wmax = 0.9;
	double wmin = 0.4;
	double dt = 1.0;
	double percent = 0.2;

	PSO pso = PSO(dim, m, Tmax, max, min, c1, c2, wmax, wmin, dt, percent);
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.updateAllParticles();
		cout << ":";
		//fitness.push_back(pso.getFitness());
		if ((1 + i) % 10 == 0) {
			fitness.push_back(pso.getFitness());

		}
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;



	}
	for (int i = 0; i < fitness.size(); i++)
	{
		cout << fitness[i] << " ";
	}
	system("pause");
}