from pathlib import Path
import pandas as pd
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel
)

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
ssa_results_dir = tinception_dir / "ssa_results"
df_lambda = pd.read_csv(ssa_results_dir / "lambdas.csv")
df_grad = pd.read_csv(ssa_results_dir / "grads.csv")

for n_component in range(1, 4):

    df_g = df_grad.query(f"component == {n_component}").copy()
    demographic_cols = ["age", "sex", "PTA", "TIV"]
    start_idx = 9
    response_cols = df_g.columns[start_idx:].tolist()

    ## run norm model
    kwargs = {
                "covariates": demographic_cols,
                "batch_effects": ["site"],
                "response_vars": response_cols, 
                "subject_ids": "subjects"
                }
    saving_dir = tinception_dir / "ssa_results" / "norm_models" / f"gradient_{n_component}"
    saving_dir.mkdir(exist_ok=True)

    ## create norm objects
    norm_train_all = NormData.from_dataframe(
                                            name="train",
                                            dataframe=df_g.query('group == "CO"'),
                                            **kwargs
                                            )
    norm_test_all = NormData.from_dataframe(
                                            name="test",
                                            dataframe=df_g,
                                            **kwargs
                                            )

    ## define models
    template_blr = BLR(
                        name="payam_blr",
                        basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
                        fixed_effect=True,
                        heteroskedastic=True,
                        warp_name="warpsinharcsinh"
                        )

    model = NormativeModel(
                    template_regression_model=template_blr,
                    savemodel=True,
                    evaluate_model=True,
                    saveresults=True,
                    saveplots=False,
                    save_dir=str(saving_dir),
                    inscaler="standardize",
                    outscaler="none",
                    )
    model.fit_predict(norm_train_all, norm_test_all)