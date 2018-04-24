//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef L2Projection_h
#define L2Projection_h

#include "MathStatement.h"
#include "DataTypes.h"

class L2Projection : public MathStatement
{
    Matrix projection;
    
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
public:
    L2Projection();
    
    L2Projection(Matrix &perm);
    
    L2Projection(const L2Projection &copy);
    
    L2Projection &operator=(const L2Projection &copy) const;
    
    virtual L2Projection *Clone();
    
    virtual ~L2Projection();
    
    Matrix GetProjectionMatrix() const;
    
    void SetProjectionMatrix(const Matrix &proj);

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
#endif /* L2Projection_h */
