#include <iostream>
#include <cmath>
#include <vector>
#include <random>

struct Option
{
	double S0;
	double K;
	double T;
	double r;
	double sigma;
	int num_simulations;
	int num_steps;
};

std::vector<double> simulate_asset_path(const Option &opt)
{
	std::vector<double> path(opt.num_steps);
	double dt = opt.T / opt.num_steps;
	double drift = (opt.r - 0.5 * opt.sigma * opt.sigma) * dt;
	double diffusion = opt.sigma * sqrt(dt);

	std::random_device rd;
	std::mt19937 generator(rd());
	std::normal_distribution<> distribution(0.0, 1.0);

	path[0] = opt.S0;
	for (int i = 1; i < opt.num_steps; ++i)
	{
		double z = distribution(generator);
		path[i] = path[i - 1] * exp(drift + diffusion * z);
	}

	return path;
}

double calculate_payoff (const Option &opt, const std::vector<double> &path)
{
	double S_T = path.back();
	return std::max(S_T - opt.K, 0.0);
}

double monte_carlo_pricing (const Option &opt)
{
	double payoff_sum = 0.0;

	for (int i = 0; i < opt.num_simulations; i++)
	{
		std::vector<double> path = simulate_asset_path(opt);
		double payoff = calculate_payoff(opt, path);
		payoff_sum += payoff;
	}

	double discount_factor = exp(-opt.r * opt.T);
	return discount_factor * (payoff_sum / opt.num_simulations);
}

int main() {

	Option opt;
    opt.S0 = 100.0;
    opt.K = 100.0;
    opt.T = 1.0;
    opt.r = 0.05;
    opt.sigma = 0.2;
    opt.num_simulations = 1000000;
    opt.num_steps = 252;

    double option_price = monte_carlo_pricing(opt);
    std::cout << "Option price: " << option_price << '\n';

	return 0;
}
