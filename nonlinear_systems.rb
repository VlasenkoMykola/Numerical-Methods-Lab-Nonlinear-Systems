#!/usr/bin/ruby -w
# -*- coding: utf-8 -*-

class NewtonSolver
  def initialize(epsilon = 1e-6, max_iterations = 100)
    @epsilon = epsilon
    @max_iterations = max_iterations
  end

  def solve_system(equations, initial_guess)
    x = initial_guess.dup
    iteration = 0

    until max_residual(equations, x) < @epsilon || iteration >= @max_iterations
      jacobian_matrix = compute_jacobian_matrix(equations, x)
      residual_vector = compute_residual_vector(equations, x)

      delta = gaussian_elimination(jacobian_matrix, residual_vector)

      x = vector_addition(x, delta)
      iteration += 1
    end

    x
  end

  private

  def compute_residual_vector(equations, x)
    #puts "in compute_residual_vector:"
    #p x
    residuals = equations.map { |eq| -eq.call(x) }
    #p residuals
  end

  def compute_jacobian_matrix(equations, x)
    num_equations = equations.size
    num_variables = x.size

    jacobian_matrix = Array.new(num_equations) { Array.new(num_variables) }

    num_equations.times do |i|
      num_variables.times do |j|
        jacobian_matrix[i][j] = derivative(equations[i], x, j)
      end
    end

    jacobian_matrix
  end

  def derivative(func, x, variable_index)
    h = 1e-8
    x_h = x.dup
    x_h[variable_index] += h

    (func.call(x_h) - func.call(x)) / h
  end

  def gaussian_elimination(matrix, vector)
    n = vector.size

    (0..n - 1).each do |i|
      max_row = i
      (i + 1..n - 1).each { |j| max_row = j if matrix[j][i].abs > matrix[max_row][i].abs }

      matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
      vector[i], vector[max_row] = vector[max_row], vector[i]

      (i + 1..n - 1).each do |j|
        factor = matrix[j][i] / matrix[i][i]

        (i..n - 1).each { |k| matrix[j][k] -= factor * matrix[i][k] }
        vector[j] -= factor * vector[i]
      end
    end

    solution = Array.new(n)

    (n - 1).downto(0) do |i|
      solution[i] = vector[i]

      (i + 1..n - 1).each { |j| solution[i] -= matrix[i][j] * solution[j] }
      solution[i] /= matrix[i][i]
    end

    solution
  end

  def vector_addition(vec1, vec2)
    vec1.zip(vec2).map { |a, b| a + b }
  end

  def max_residual(equations, x)
    residuals = compute_residual_vector(equations, x)
    #puts "in max_residual:"
    #p residuals
    residuals.map(&:abs).max
  end
end

solver = NewtonSolver.new

# System #1 of non-linear equations
system1 = [
  ->(x) { x[0] * x[0] - 2 * x[0] * x[1] + 1 },
  ->(x) { x[0] * x[0] + x[1] * x[1] - 2 }
]

#p system1

initial_guess = [-0.5, -0.5]
solution = solver.solve_system(system1, initial_guess)
puts "System 1: initial: #{initial_guess} solution: #{solution}"

initial_guess = [1.1, 1.1]
solution = solver.solve_system(system1, initial_guess)
puts "System 1: initial: #{initial_guess} solution: #{solution}"

# System #2 of non-linear equations
n = 4
system2 = Array.new(n) { |i| ->(x) {
    sum = 0.0 - n
    for j in 0..n-1 do
      if i = j then
        sum = sum + x[j] * x[j] * x[j]
      else
        sum = sum + x[j] * x[j]
      end
    end
    sum
  }
}

#p system2


#initial_guess = Array.new(n, 1.0)
#solution = solver.solve_system(system2, initial_guess)
#puts "System 2: initial: #{initial_guess} solution: #{solution}"

initial_guess = Array.new(n, 0.99999999)
solution = solver.solve_system(system2, initial_guess)
puts "System 2: initial: #{initial_guess} solution: #{solution}"
