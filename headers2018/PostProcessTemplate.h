//
//  PostProcessTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//

#ifndef PostProcessTemplate_h
#define PostProcessTemplate_h

//#include "PostProcess.h"
#include <list>

template<class math>
class PostProcessTemplate
{
    std::vector<typename math::PostProcVar> scalarvariables;
    std::vector<typename math::PostProcVar> vectorvariables;
    
    AppendVariable(typename math::postprocvar);
    
    std::vector<double> PostProcResult(MathStatement &mathStatement, unsigned int varIndex, IntPontData &data) const {
        math *locptr = dynamic_cast<typename math*> &mathStatement;
        if(!locptr) DebugStop();
        const numScalarVariables = NumScalarVariables();
        return locptr->PostProcessSolution(varIndex < numScalarVariables? scalarvariable[varIndex] : vectorvariables[varIndex-numScalarVariables], data);
    }
    
    inline unsigned int NumScalarVariables() const {
        return scalarvariables.size();
    }
    
    inline unsigned int NumVectorVariables() const {
        return vectorvariables.size();
    }
    
    inline unsigned int NumVariables() const {
        return NumScalarVariables() + NumVectorVariables();
    }
};
#endif /* PostProcessTemplate_h */
