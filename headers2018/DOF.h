//
//  DOF.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef DOF_h
#define DOF_h

class DOF
{
    int64_t firstequation = -1;
    
    int nshape = 0;
    
    int nstate = 0;
    
public:
    
    DOF();
    
    DOF(const DOF &copy);
    
    DOF &operator=(const DOF &copy);
    
    ~DOF();
    
    int64_t GetFirstEquation();
    
    void SetFirstEquation(int64_t first);
    
    void SetNShapeState(int NShape, int NState);
    
    int GetNShape() const;
    
    int GetNState() const;
    
    
};
#endif /* DOF_h */
