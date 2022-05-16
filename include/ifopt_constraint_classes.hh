
#include "includes_and_types.hh"

#pragma once

namespace ifopt {

class NullSpaceRepresentationConstraint : public ConstraintSet {
    // MEPPF 94
    NullSpaceRepresentationConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class SteadyStateConstraint : public ConstraintSet {
    // MEPPF 95
public:
    SteadyStateConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class SmoothConstraint : public ConstraintSet {
    // MEPPF 96
public:
    SmoothConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class RelaxedFluxUpperConstraint : public ConstraintSet {
    // MEPPF 97
public:
    RelaxedFluxUpperConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class RelaxedFluxLowerConstraint : public ConstraintSet {
    // MEPPF 98
public:
    RelaxedFluxLowerConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class SignConstraint : public ConstraintSet {
    // MEPPF 99
public:
    SignConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class RelaxedFluxSignConstraint : public ConstraintSet {
    // MEPPF 100
public:
    RelaxedFluxSignConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class MetabolitesUpperBoundConstraint : public ConstraintSet {
    // MEPPF 101 upper
public:
    MetabolitesUpperBoundConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

class MetabolitesLowerBoundConstraint : public ConstraintSet {
    // MEPPF 101 lower
public:
    MetabolitesLowerBoundConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
}

} // namespace ifopt