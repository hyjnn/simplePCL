#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <format>
#include <numbers>

#include "FDTD.hpp"
#include "MDVector.hpp"

void fundamentalMode(PIC::FieldSolver &sim) { // Assumes square sim
    auto step = sim.getSpaceStep();
    auto max_index = sim.getShape()[0];
    auto time_step = step * sim.getStepRatio();
    auto length = step * (max_index - 1);
    auto frequency = std::numbers::pi * PIC::c / (length);
    auto& fields = sim.getFields();

    PIC::floatType phase, Ey0, Hz0;
    for (std::size_t i = 0; i < max_index; i++) {
        phase = std::numbers::pi * i * step/length;
        Ey0 = -PIC::c * PIC::mu0 * std::sin(phase); // Hz0 = std::cos(phase) * std::sin(frequency*time_step);
        for (std::size_t j = 0; j < max_index; j++) {
            fields[1][i, j] = Ey0;
        }
    }
}



int main() {
    PIC::SimEngine sim({ 200, 200 }, 2, 1e-4);
    auto &field_sim = sim.getFieldSim();
    auto &particle_sim = sim.getParticleSim();
    auto &positions = particle_sim.getPositions();
    fundamentalMode(field_sim);

    positions[0, 0] = 100.5 * field_sim.getSpaceStep();
    positions[0, 1] = 101.5 * field_sim.getSpaceStep();
    positions[1, 0] = 100.5 * field_sim.getSpaceStep();
    positions[1, 1] = 98.5 * field_sim.getSpaceStep();
    particle_sim.getMasses()[0] = 9.109e-31;
    particle_sim.getCharges()[0] = 1.602e-13;
    particle_sim.getMasses()[1] = 9.109e-31;
    particle_sim.getCharges()[1] = 1.602e-13;

    sim.initialize();

    field_sim.exportToFile(std::format("field_data{}.txt", 0));
    for (int i = 1; i < 10; i++) {
        sim.run(50);
        std::cout << std::format("Batch {} done\n", i);
        field_sim.exportToFile(std::format("field_data{}.txt", i));
    }
    sim.exportTracked("particle_data.txt");
}