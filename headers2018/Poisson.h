//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef Poisson_h
#define Poisson_h

#include "MathStatement.h"
#include "DataTypes.h"

class Poisson : public MathStatement
{
    Matrix permeability;
    
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
public:
    
    Poisson();
    
    Poisson(Matrix &perm);
    
    Poisson(const Poisson &copy);
    
    Poisson &operator=(const Poisson &copy) const;
    
    virtual Poisson *Clone();

    virtual ~Poisson();
    
    Matrix GetPermeability() const;
    
    void SetPermeability(const Matrix &perm);
    
    std::function<void(const VecDouble &co, VecDouble &result)> GetForceFunction() const
    {
        return forceFunction;
    }
    
    void SetForceFunction(const std::function<void(const VecDouble &co, VecDouble &result)> &f)
    {
        forceFunction = f;
    }
    
    virtual int NState() const;
    
    virtual void Contribute(IntPointData &integrationpointdata, Matrix &EK, Matrix &EF) const;
    
    virtual void PostProcess(IntPointData &integrationpointdata, const std::string &variable,
                             VecDouble &postprocvalue) const;

};
#endif /* Poisson_h */
