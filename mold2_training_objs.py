from dataclasses import dataclass

@dataclass
class BaseClassifierTrainingObjects:
  classifier_name: str # "knn", ..., "xgboost"
  classifiers: list[list[any]]
  seeds: list[list[str]]
  scalers:list[list[any]]

@dataclass
class TrainingObjects:
  features: list[str]
  best_model: any
  mcc: any
  scaler: any
  base_cl_objs: dict[str, BaseClassifierTrainingObjects]
