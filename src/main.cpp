#include <cmath>
#include <iostream>
#include <random>
#include <vector>

class Particle {
  int index;
  double position;
  double value;
  double velocity;
  double c1, c2;

  double _best_pos;
  double _best_value;

  double (*f)(double);

public:
  Particle(int _index, double _position, double _velocity, double _c1,
           double _c2, double (*_f)(double))
      : index(_index), position(_position), velocity(_velocity), c1(_c1),
        c2(_c2), f(_f) {
    value = f(position);
    _best_pos = position;
    _best_value = f(_best_pos);
  }

  void step(double step_weight, double r1, double r2, double best_global_pos) {
    velocity = step_weight * velocity + c1 * r1 * (_best_pos - position) +
               c2 * r2 * (best_global_pos - position);
    position = position + velocity;

    value = f(position);

    if (value < _best_value) {
      _best_value = value;
      _best_pos = position;
    }
  }

  double best_value() { return _best_value; }

  double best_pos() { return _best_pos; }

  double get_position() { return position; }

  double get_value() { return value; }

  int get_index() { return index; }

  void print() {
    std::cout << "Particle " << index << std::endl
              << "Position: " << position << std::endl
              << "Velocity: " << velocity << std::endl
              << "Best Pos: " << _best_pos << std::endl
              << "Best Value: " << _best_value << std::endl;
  }
};

class PSO {
  double best_global_pos;
  double best_global_value;
  int steps;

  double c1, c2;
  double (*f)(double);
  int n_particles;

  std::vector<Particle> particles;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform{0.0, 1.0};

public:
  PSO(double initial_pos, double _c1, double _c2, double (*_f)(double),
      int _n_particles)
      : c1(_c1), c2(_c2), f(_f), n_particles(_n_particles) {
    steps = 0;

    best_global_pos = initial_pos;
    best_global_value = f(best_global_pos);
    particles.reserve(n_particles);

    double position, velocity;
    for (int i = 0; i < n_particles; ++i) {
      position = uniform(generator);
      velocity = uniform(generator);
      particles.push_back(Particle(i, position, velocity, c1, c2, f));
    }
  }

  void step(int n) {
    double step_weight = 1.0 / (n + 1);

    double r1 = uniform(generator);
    double r2 = uniform(generator);
    for (auto &p : particles) {
      p.step(step_weight, r1, r2, best_global_pos);

      if (p.best_value() < best_global_value) {
        best_global_pos = p.best_pos();
        best_global_value = p.best_value();
      }
    }
  }

  void print_global() {
    std::cout << "Particle Max Dist: " << particle_max_dist() << std::endl
              << "Particle Value Max Dist: " << particle_value_max_dist()
              << std::endl
              << "Steps: " << steps << std::endl
              << "Global Best Pos: " << best_global_pos << std::endl
              << "Global Best Value: " << best_global_value << std::endl;
  }

  void print(int i = -1) {
    for (auto &p : particles) {
      if (i > 0) {
        if (p.get_index() == i) {
          p.print();
          print_global();
        }
      } else {
        p.print();
        print_global();
      }
    }
  }

  double particle_max_dist() {
    double min_pos = particles[0].get_position();
    double max_pos = particles[0].get_position();

    for (auto &p : particles) {
      if (p.get_position() < min_pos)
        min_pos = p.get_position();
      else if (p.get_position() > max_pos)
        max_pos = p.get_position();
    }
    return max_pos - min_pos;
  }

  double particle_value_max_dist() {
    double min_value = particles[0].get_value();
    double max_value = particles[0].get_value();

    for (auto &p : particles) {
      if (p.get_value() < min_value)
        min_value = p.get_value();
      else if (p.get_value() > max_value)
        max_value = p.get_value();
    }
    return max_value - min_value;
  }

  void run(double threshold, int max_steps) {
    while (steps < max_steps) {
      step(steps);
      ++steps;
      if ((particle_value_max_dist() < threshold) &&
          (particle_max_dist() < threshold))
        break;
    }
  }
};

double f(double x) { return pow(x, 2) - 3 * x + 1; }

int main() {
  PSO pso(1.0, 0.5, 0.5, f, 100);
  pso.run(0.0001, 50);
  pso.print_global();
}
