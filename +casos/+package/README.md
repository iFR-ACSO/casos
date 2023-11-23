# Internal class diagram

```mermaid
flowchart LR
  subgraph casadi
    casadi.Callback[Callback] --> casadi.Function[Function]
  end
  subgraph casos
    casos.Function[Function]
    casos.conic[conic]
    casos.sdpsol[sdpsol]
    casos.sossol[sossol]
  end
  subgraph casos.package.functions
    FunctionInterface --> FunctionCommon
    FunctionWrapper -.-> FunctionInterface
    CasadiFunction --> FunctionInterface
    PSFunction --> FunctionInterface
  end
  subgraph casos.package.solvers
    ConicSolver --> SolverCallback
    SedumiInterface --> ConicSolver
    SCSInterface --> ConicSolver
    SdpsolInternal --> SolverCallback
    SdpsolInternal -.-> ConicSolver
    conicInternal -.-> SedumiInterface
    conicInternal -.-> SCSInterface
    SossolInternal --> FunctionWrapper
    SossolInternal -.-> SossdpRelaxation
    SossdpRelaxation --> FunctionInterface
    SossdpRelaxation -.-> SdpsolInternal
  end
  casos.Function --> FunctionWrapper
  casos.Function -.-> CasadiFunction
  casos.Function -.-> PSFunction
  CasadiFunction -.-> casadi.Function
  SolverCallback --> casadi.Callback
  SolverCallback --> FunctionCommon
  casos.conic -.-> conicInternal
  casos.sdpsol -.-> SdpsolInternal
  casos.sossol -.-> SossolInternal
```
