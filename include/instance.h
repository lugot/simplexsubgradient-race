#ifndef INCLUDE_INSTANCE_H_
#define INCLUDE_INSTANCE_H_

#include <ostream>

class Instance {
   public:
    Instance();
    // wrapper of IloAlgorithm::Status
    enum class SolutionStatus {
        Unknown,
        Feasible,
        Optimal,
        Infeasible,
        Unbounded,
        InfeasibleOrUnbounded,
        Error,
        EmptyInstance,
        NeverSolved,
        ModelChanged,
    };

    // TODO(lugot): ADD virtual thinghy for solve
    SolutionStatus getStatus();
    bool isSolved();

   protected:
    SolutionStatus status;
};

std::ostream& operator<<(std::ostream& os, const Instance::SolutionStatus& s);

#endif  // INCLUDE_INSTANCE_H_
