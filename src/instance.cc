#include "../include/instance.h"

#include "ilcplex/ilocplexi.h"

Instance::Instance() { status = SolutionStatus::NeverSolved; }

Instance::SolutionStatus Instance::getStatus() { return status; }

bool Instance::isSolved() {
    return status == SolutionStatus::Feasible ||
           status == SolutionStatus::Optimal ||
           status == SolutionStatus::Infeasible ||
           status == SolutionStatus::Unbounded ||
           status == SolutionStatus::InfeasibleOrUnbounded;
}

std::ostream& operator<<(std::ostream& os, const Instance::SolutionStatus& s) {
    switch (s) {
        case Instance::SolutionStatus::EmptyInstance:
            os << "EmptyInstance";
            break;
        case Instance::SolutionStatus::NeverSolved:
            os << "NeverSolved";
            break;
        case Instance::SolutionStatus::ModelChanged:
            os << "ModelChanged";
            break;
        default:
            // get back to IloAlgorithm::Status << overloaded operator
            os << static_cast<IloAlgorithm::Status>(s);
            break;
    }

    return os;
}
