/*
 An optimization problem.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_ASLPROBLEM_H_
#define MP_ASLPROBLEM_H_

#include <algorithm>
#include <deque>
#include <vector>

#include "aslexpr.h"

namespace mp {

template <typename SuffixPtr>
class SuffixData;

namespace asl {
namespace internal {
class ASLBuilder;
}
}

// A solution of an optimization problem.
class Solution {
 private:
  int solve_code_;
  int num_vars_;
  int num_cons_;
  double *values_;
  double *dual_values_;

  FMT_DISALLOW_COPY_AND_ASSIGN(Solution);

 public:
  // Constructs a solution with zero variables and constraints and the
  // solve code -1.
  Solution();

  ~Solution();

  // Swaps this solution with other.
  void Swap(Solution &other);

  // Returns the solution status.
  sol::Status status() const {
    return solve_code_ < 0 || solve_code_ >= 600 ? sol::UNKNOWN
        : static_cast<sol::Status>(solve_code_ / 100 * 100);
  }

  // Returns the solve code.
  int solve_code() const { return solve_code_; }

  // Returns the number of variables.
  int num_vars() const { return num_vars_; }

  // Returns the number of constraints.
  int num_cons() const { return num_cons_; }

  // Returns the values of all variables.
  const double *values() const { return values_; }

  // Returns the values of all dual variables.
  const double *dual_values() const { return dual_values_; }

  // Returns the value of a variable.
  double value(int var) const {
    assert(var >= 0 && var < num_vars_);
    return values_[var];
  }

  // Returns the value of a dual variable corresponding to constraint con.
  double dual_value(int con) const {
    assert(con >= 0 && con < num_cons_);
    return dual_values_[con];
  }

  // Reads a solution from the file <stub>.sol.
  void Read(fmt::CStringRef stub, int num_vars, int num_cons);
};

class ASLSuffixPtr {
 private:
  class Proxy {
   private:
    ASL *asl_;
    SufDesc *ptr_;

    friend class ASLSuffixPtr;
    friend class SuffixView;

    template <typename T, typename Visitor>
    static void VisitValues(
        int num_items, const T *values, const int *map, Visitor &visitor);

    Proxy(ASL *asl, SufDesc *p) : asl_(asl), ptr_(p) {}

    friend class mp::SuffixData<ASLSuffixPtr>;

    void set_data(int *values, std::size_t) {
      ptr_->u.i = values;
    }

   public:
    const char *name() const { return ptr_->sufname; }
    int kind() const { return ptr_->kind; }

    bool has_values() const { return ptr_->u.i != 0; }

    int int_value(int index) const { return ptr_->u.i[index]; }
    void set_value(int index, int value) { ptr_->u.i[index] = value; }

    // Iterates over nonzero suffix values and sends them to the visitor.
    template <typename Visitor>
    void VisitValues(Visitor &visitor) const;
  };

  mutable Proxy proxy_;

  friend class ASLProblem;
  friend class SuffixView;

  ASLSuffixPtr(ASL *asl, SufDesc *s) : proxy_(asl, s) {}

  void True() const {}
  typedef void (ASLSuffixPtr::*SafeBool)() const;

 public:
  // Constructs an ASLSuffixPtr object representing a null pointer.
  ASLSuffixPtr() : proxy_(0, 0) {}

  operator SafeBool() const {
    return proxy_.ptr_ ? &ASLSuffixPtr::True : 0;
  }

  Proxy *operator->() const { return &proxy_; }

  Proxy operator*() const { return proxy_; }

  // TODO
  void set_value(int index, int value) { proxy_.set_value(index, value); }
};

template <typename SuffixType>
SuffixType Cast(ASLSuffixPtr s) { return s; } // TODO

template <typename T, typename Visitor>
void ASLSuffixPtr::Proxy::VisitValues(
    int num_items, const T *values, const int *map, Visitor &visitor) {
  if (!values) return;
  if (map) {
    for (int i = 0; i < num_items; ++i) {
      int index = map[i];
      T value = values[i];
      if (value != 0 && index >= 0)
        visitor.Visit(index, value);
    }
  } else {
    for (int i = 0; i < num_items; ++i) {
      if (T value = values[i])
        visitor.Visit(i, value);
    }
  }
}

template <typename Visitor>
void ASLSuffixPtr::Proxy::VisitValues(Visitor &visitor) const {
  int kind = ptr_->kind & mp::internal::SUFFIX_KIND_MASK;
  int num_items = (&asl_->i.n_var_)[kind];
  const int *map = kind < 2 ? (&asl_->i.vmap)[kind] : 0;
  if ((ptr_->kind & suf::FLOAT) != 0)
    VisitValues(num_items, ptr_->u.r, map, visitor);
  else
    VisitValues(num_items, ptr_->u.i, map, visitor);
}

// A view on a collection of suffixes.
// Unlike a real collection, a view doesn't own elements, but only provides
// access to them, like an iterator. In fact, the view can be considered as
// a pair of iterators. Therefore a SuffixView object should only be used
// while the ASLProblem object containing suffixes is alive and hasn't been
// modified since the view was obtained.
class SuffixView {
 private:
  ASL *asl_;
  int kind_;

  friend class ASLProblem;
  friend class asl::internal::ASLBuilder;

  SuffixView(ASL *asl, int kind) : asl_(asl), kind_(kind) {}

 public:
  SuffixView() : asl_(0), kind_(0) {}

  class iterator :
    public std::iterator<std::forward_iterator_tag, ASLSuffixPtr::Proxy> {
   private:
    ASL *asl_;
    SufDesc *ptr_;

   public:
    iterator() : asl_(0), ptr_(0) {}
    iterator(ASL *asl, SufDesc *p) : asl_(asl), ptr_(p) {}

    ASLSuffixPtr::Proxy operator*() const {
      return ASLSuffixPtr::Proxy(asl_, ptr_);
    }

    ASLSuffixPtr operator->() const { return ASLSuffixPtr(asl_, ptr_); }

    iterator &operator++() {
      ptr_ = ptr_->next;
      return *this;
    }

    iterator operator++(int ) {
      iterator it(*this);
      ++*this;
      return it;
    }

    bool operator==(iterator other) const { return ptr_ == other.ptr_; }
    bool operator!=(iterator other) const { return ptr_ != other.ptr_; }
  };

  // Finds a suffix with specified name.
  ASLSuffixPtr Find(const char *name, unsigned flags = 0) const;

  iterator begin() const { return iterator(asl_, asl_->i.suffixes[kind_]); }
  iterator end() const { return iterator(); }
};

class ProblemChanges;

// An optimization problem.
class ASLProblem {
 private:
  ASL *asl_;
  int var_capacity_;
  int obj_capacity_;
  int logical_con_capacity_;

  // Array of variable types or null if continuous variables precede
  // integer and binary variables.
  var::Type *var_types_;

  FMT_DISALLOW_COPY_AND_ASSIGN(ASLProblem);

  static void IncreaseCapacity(int size, int &capacity) {
    if (capacity == 0 && size != 0)
      throw Error("Problem can't be modified");
    capacity = (std::max)(capacity, size);
    capacity = capacity ? 2 * capacity : 8;
  }

  template <typename T>
  static void Grow(T *&array, int &size, int &capacity) {
    T *new_array = new T[capacity];
    std::copy(array, array + size,
      fmt::internal::make_ptr(new_array, capacity));
    delete [] array;
    array = new_array;
  }

  friend class asl::internal::ASLBuilder;
  friend class ASLSolver;

  // Frees all the arrays that were allocated by modifications to the problem.
  void Free();

  // Write an .nl file.
  void WriteNL(fmt::CStringRef stub, ProblemChanges *pc = 0,
               unsigned flags = 0);

  class Proxy {
   private:
    mutable ASL *asl_;
    unsigned flags_;

    Proxy(ASL *asl, unsigned flags = 0) : asl_(asl), flags_(flags) {}

    friend class ASLProblem;
    friend class ASLSolver;
    friend class asl::internal::ASLBuilder;

    void Free() {
      if (asl_)
        ASL_free(&asl_);
    }

   public:
    Proxy(const Proxy &other) : asl_(other.asl_), flags_(other.flags_) {
      other.asl_ = 0;
    }
    Proxy &operator=(const Proxy &other) {
      Free();
      asl_ = other.asl_;
      other.asl_ = 0;
      flags_ = other.flags_;
      return *this;
    }
    ~Proxy() { Free(); }
  };

  static double lb(int index, int size, const double *lbs, const double *ubs) {
    internal::Unused(size);
    assert(index >= 0 && index < size);
    return lbs[ubs ? index : (index * 2)];
  }

  static double ub(int index, int size, const double *lbs, const double *ubs) {
    internal::Unused(size);
    assert(index >= 0 && index < size);
    return ubs ? ubs[index] : lbs[index * 2 + 1];
  }

  template <typename ProblemType>
  struct BasicProblemItem {
    ProblemType *problem_;
    int index_;

    typedef ProblemType Problem;

    BasicProblemItem(ProblemType *p, int index)
      : problem_(p), index_(index) {}
  };

  typedef BasicProblemItem<const ASLProblem> ProblemItem;
  typedef BasicProblemItem<ASLProblem> MutProblemItem;

 public:
  typedef asl::internal::ExprTypes ExprTypes;

  ASLProblem();
  ASLProblem(Proxy proxy);

  ~ASLProblem();

  const char *name() const { return asl_->i.filename_; }

  // Returns the number of variables.
  int num_vars() const { return asl_->i.n_var_; }

  // Returns the number of objectives.
  int num_objs() const { return asl_->i.n_obj_; }

  // Returns the number of algebraic constraints.
  int num_algebraic_cons() const { return asl_->i.n_con_; }

  // Returns the number of integer variables including binary.
  int num_integer_vars() const {
    return asl_->i.nbv_ + asl_->i.niv_ + asl_->i.nlvbi_ +
        asl_->i.nlvci_ + asl_->i.nlvoi_;
  }

  // Returns the number of continuous variables.
  int num_continuous_vars() const {
    return num_vars() - num_integer_vars();
  }

  // Returns the number of nonlinear objectives.
  int num_nonlinear_objs() const { return asl_->i.nlo_; }

  // Returns the number of nonlinear constraints.
  int num_nonlinear_cons() const { return asl_->i.nlc_; }

  // Returns the number of logical constraints.
  int num_logical_cons() const { return asl_->i.n_lcon_; }

  // Returns the number of nonzeros in constraints' Jacobian.
  int num_con_nonzeros() const { return asl_->i.nzc_; }

  // Returns the number of common expressions.
  int num_common_exprs() const { return asl_->i.ncom0_ + asl_->i.ncom1_; }

  // Returns the type of the variable.
  var::Type var_type(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    if (var_types_)
      return var_types_[var_index];
    return var_index < num_continuous_vars() ? var::CONTINUOUS : var::INTEGER;
  }

  // An optimization variable.
  class Variable : private ProblemItem {
   private:
    friend class ASLProblem;

    Variable(const ASLProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const ASLProblem &p) {
      return p.num_vars();
    }

   public:
    // Returns the index of the variable.
    int index() { return index_; }

    // Returns the lower bound on the variable.
    double lb() const {
      return problem_->asl_->i.LUv_[
          problem_->asl_->i.Uvx_ ? index_ : (index_ * 2)];
    }

    // Returns the upper bound on the variable.
    double ub() const {
      return problem_->asl_->i.Uvx_ ?
            problem_->asl_->i.Uvx_[index_] :
            problem_->asl_->i.LUv_[index_ * 2 + 1];
    }

    // Returns the type of the variable.
    var::Type type() const {
      if (problem_->var_types_)
        return problem_->var_types_[index_];
      return index_ < problem_->num_continuous_vars() ?
            var::CONTINUOUS : var::INTEGER;
    }

    bool operator==(Variable other) const {
      MP_ASSERT(problem_ == other.problem_,
                "comparing variables from different problems");
      return index_ == other.index_;
    }
    bool operator!=(Variable other) const {
      return !(*this == other);
    }
  };

  // Returns the variable at the specified index.
  Variable var(int index) const {
    MP_ASSERT(0 <= index && index < num_vars(), "invalid index");
    return Variable(this, index);
  }

  // Returns the lower bounds for the variables.
  const double *var_lb() const { return asl_->i.LUv_; }

  // Returns the upper bounds for the variables.
  const double *var_ub() const { return asl_->i.Uvx_; }

  // Returns the initial values for the variables.
  const double *initial_values() const { return asl_->i.X0_; }

  // Returns the lower bounds for the constraints.
  const double *con_lb() const { return asl_->i.LUrhs_; }

  // Returns the upper bounds for the constraints.
  const double *con_ub() const { return asl_->i.Urhsx_; }

  // An objective.
  class Objective : private ProblemItem {
   private:
    friend class ASLProblem;

    Objective(const ASLProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const ASLProblem &p) {
      return p.num_objs();
    }

   public:
    // Returns the type of the objective.
    obj::Type type() const {
      return static_cast<obj::Type>(problem_->asl_->i.objtype_[index_]);
    }

    // Returns the linear part of the objective expression.
    asl::LinearObjExpr linear_expr() const {
      return asl::LinearObjExpr(problem_->asl_->i.Ograd_[index_]);
    }

    // Returns the nonlinear part of the objective expression.
    asl::NumericExpr nonlinear_expr() const {
      if (problem_->asl_->i.ASLtype != ASL_read_fg)
        return asl::NumericExpr();
      return asl::NumericExpr::Create<asl::NumericExpr>(
            reinterpret_cast<ASL_fg*>(problem_->asl_)->I.obj_de_[index_].e);
    }

    bool operator==(Objective other) const {
      MP_ASSERT(problem_ == other.problem_,
                "comparing objectives from different problems");
      return index_ == other.index_;
    }

    bool operator!=(Objective other) const {
      return !(*this == other);
    }
  };

  // Returns the objective at the specified index.
  Objective obj(int index) const {
    MP_ASSERT(0 <= index && index < num_objs(), "invalid index");
    return Objective(this, index);
  }

  // An algebraic constraint.
  class AlgebraicCon : private ProblemItem {
   private:
    friend class ASLProblem;

    AlgebraicCon(const ASLProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const ASLProblem &p) {
      return p.num_objs();
    }

   public:
    // Returns the lower bound on the constraint.
    double lb() const {
      return problem_->asl_->i.LUrhs_[
          problem_->asl_->i.Urhsx_ ? index_ : (index_ * 2)];
    }

    // Returns the upper bound on the constraint.
    double ub() const {
      return problem_->asl_->i.Urhsx_ ?
            problem_->asl_->i.Urhsx_[index_] :
            problem_->asl_->i.LUrhs_[index_ * 2 + 1];
    }

    // Returns the linear part of the constraint expression.
    asl::LinearConExpr linear_expr() const {
      return asl::LinearConExpr(problem_->asl_->i.Cgrad_[index_]);
    }

    // Returns the nonlinear part of the constraint expression.
    asl::NumericExpr nonlinear_expr() const {
      if (problem_->asl_->i.ASLtype != ASL_read_fg)
        return asl::NumericExpr();
      return asl::NumericExpr::Create<asl::NumericExpr>(
            reinterpret_cast<ASL_fg*>(problem_->asl_)->I.con_de_[index_].e);
    }

    bool operator==(AlgebraicCon other) const {
      MP_ASSERT(problem_ == other.problem_,
                "comparing constraint from different problems");
      return index_ == other.index_;
    }

    bool operator!=(AlgebraicCon other) const {
      return !(*this == other);
    }
  };

  // Returns the algebraic constraint at the specified index.
  AlgebraicCon algebraic_con(int index) const {
    MP_ASSERT(0 <= index && index < num_algebraic_cons(), "invalid index");
    return AlgebraicCon(this, index);
  }

  class CommonExpr : private ProblemItem {
   private:
    friend class ASLProblem;

    CommonExpr(const ASLProblem *p, int index) : ProblemItem(p, index) {}

    bool IsSupported() const {
      return problem_->asl_->i.ASLtype == ASL_read_fg;
    }

   public:
    // Returns the linear part of the common expression.
    asl::LinearCommonExpr linear_expr() const {
      if (!IsSupported()) return asl::LinearCommonExpr();
      return asl::LinearCommonExpr(
            &reinterpret_cast<ASL_fg*>(problem_->asl_)->I.cexps_[index_]);
    }

    // Returns the nonlinear part of the common expression.
    asl::NumericExpr nonlinear_expr() const {
      if (!IsSupported()) return asl::NumericExpr();
      return asl::NumericExpr::Create<asl::NumericExpr>(
            reinterpret_cast<ASL_fg*>(problem_->asl_)->I.cexps_[index_].e);
    }
  };

  // Returns the common expression at the specified index.
  CommonExpr common_expr(int index) const {
    return CommonExpr(this, index);
  }

  template <typename ExprType>
  ExprType GetExpr(cde *Edag1info::*ptr, int index, int size) const {
    internal::Unused(size);
    assert(index >= 0 && index < size);
    if (asl_->i.ASLtype != ASL_read_fg)
      return ExprType();
    return asl::Expr::Create<ExprType>(
          (reinterpret_cast<ASL_fg*>(asl_)->I.*ptr)[index].e);
  }
  asl::NumericExpr GetExpr(cde *Edag1info::*ptr, int index, int size) const {
    return GetExpr<asl::NumericExpr>(ptr, index, size);
  }

  // Returns a logical constraint expression.
  asl::LogicalExpr logical_con_expr(int lcon_index) const {
    return GetExpr<asl::LogicalExpr>(
          &Edag1info::lcon_de_, lcon_index, num_logical_cons());
  }

  // Returns the name of the variable or a null pointer if the name is not
  // available.
  const char *var_name(int var_index) const {
    return var_name_ASL(asl_, var_index);
  }

  // Returns the name of the constraint or a null pointer if the name is not
  // available.
  const char *con_name(int con_index) const {
    return con_name_ASL(asl_, con_index);
  }

  // Returns the name of the logical constraint or a null pointer if the
  // name is not available.
  const char *logical_con_name(int lcon_index) const {
    return lcon_name_ASL(asl_, lcon_index);
  }

  class ColMatrix {
   private:
    Edaginfo *info_;

   public:
    explicit ColMatrix(Edaginfo *info) : info_(info) {}

    // Returns an offset of a column.
    int col_start(int col_index) const {
      return info_->A_colstarts_[col_index];
    }

    const int *col_starts() const { return info_->A_colstarts_; }

    // Returns the row index of an element.
    int row_index(int elt_index) const { return info_->A_rownos_[elt_index]; }

    const int *row_indices() const { return info_->A_rownos_; }

    // Returns the value of an element.
    double value(int elt_index) const { return info_->A_vals_[elt_index]; }

    const double *values() const { return info_->A_vals_; }
  };

  // Returns the columnwise representation of the constraint matrix.
  ColMatrix col_matrix() const { return ColMatrix(&asl_->i); }

  // Returns the solve code.
  int solve_code() const { return asl_->p.solve_code_; }

  // Sets the solve code.
  void set_solve_code(int value) {
    asl_->p.solve_code_ = value;
  }

  SuffixView suffixes(int kind) const {
    assert(kind <= internal::NUM_SUFFIX_KINDS);
    return SuffixView(asl_, kind);
  }

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type = var::CONTINUOUS);

  // Adds an objective.
  void AddObj(obj::Type type, asl::NumericExpr expr);

  // Adds a logical constraint.
  void AddCon(asl::LogicalExpr expr);

  // Flags for the Read method.
  enum { READ_INITIAL_VALUES = 1 };

  // Reads a problem from the file <stub>.nl.
  void Read(fmt::StringRef stub, unsigned flags = 0);

  // Flags for the Solve method.
  enum { IGNORE_FUNCTIONS = 1 };

  // Solves the current problem.
  void Solve(fmt::StringRef solver_name, Solution &sol,
      ProblemChanges *pc = 0, unsigned flags = 0);
};

// Writes the linear part of the problem in the AMPL format.
fmt::Writer &operator<<(fmt::Writer &w, const ASLProblem &p);

// Changes (additions) to an optimization problem.
class ProblemChanges {
 private:
  const ASLProblem *problem_;
  std::vector<double> var_lb_;
  std::vector<double> var_ub_;
  std::vector<double> con_lb_;
  std::vector<double> con_ub_;
  std::deque<ograd> con_terms_;
  std::deque<ograd> obj_terms_;
  std::vector<ograd*> cons_;
  std::vector<ograd*> objs_;
  std::vector<char> obj_types_;
  NewVCO vco_;

  friend class ASLProblem;

  NewVCO *vco();

 public:
  explicit ProblemChanges(const ASLProblem &p) : problem_(&p), vco_() {}

  ProblemChanges(const ProblemChanges &other);
  ProblemChanges &operator=(const ProblemChanges &rhs);

  // Returns the number of additional variables.
  int num_vars() const { return static_cast<int>(var_lb_.size()); }

  // Returns the number of additional constraints.
  int num_cons() const { return static_cast<int>(cons_.size()); }

  // Returns the number of additional objectives.
  int num_objs() const { return static_cast<int>(objs_.size()); }

  // Adds a variable.
  int AddVar(double lb, double ub) {
    var_lb_.push_back(lb);
    var_ub_.push_back(ub);
    return static_cast<int>(problem_->num_vars() + var_lb_.size() - 1);
  }

  // Adds an objective.
  void AddObj(obj::Type type,
      unsigned size, const double *coefs, const int *vars);

  // Adds a constraint.
  void AddCon(const double *coefs, double lb, double ub);
};
}  // namespace mp

#endif  // MP_ASLPROBLEM_H_
