//
//  MathStatement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef MathStatement_h
#define MathStatement_h

class MathStatement
{
    int matid = 0;
    
    int nstate = 0;
    
public:
    
    MathStatement();
    
    MathStatement(const MathStatement &copy);
    
    MathStatement &operator=(const MathStatement &copy) const;
    
    virtual ~MathStatement();
    
    virtual MathStatement *Clone() = 0;
    
    virtual int NState() const = 0;
    
    virtual void Contribute(IntPointData &integrationpointdata, Matrix &EK, Matrix &EF) const = 0;
    
    virtual void PostProcess(IntPointData &integrationpointdata, const std::string &variable,
                             VecDouble &postprocvalue) const = 0;
    
};
#endif /* MathStatement_h */
