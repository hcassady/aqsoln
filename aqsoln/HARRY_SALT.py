"""Calculate some harry salt properties."""
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# Train model
# Load data and map salt identities onto number
salt_to_id = {"None": 0, "HCl": 1, "NaCl": 2, "NaOH": 3}
training_data = pd.read_csv("physical_data.csv")
training_data["salt_id"] = training_data["salt"].apply(
    lambda s: salt_to_id[s]
)

# Train data
X = training_data[["salt_id", "weight_percent", "temperature"]].values
y = training_data["density"]

# Split into test/train
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=0
)

# Test model
model = RandomForestRegressor(random_state=0).fit(X_train, y_train)
model.score(X_test, y_test)


def get_density(salt: str="None", weight_percent: float=0, temperature: float=0) -> float:
    """Returns the density of a solution of salt

    Calculates the density of a solution of salt at a specific weight percent and temperature. 

    Accepts
    -------
    salt : float
        Species of salt: None, HCl, NaCl, NaOH, by default None
    weight_percent : float
        The weight percent of the salt in solution between 0 and 1, by default 0
    temperature : float
        The temperature of the solution in Celcius, by default 0

    Returns
    -------
    float
        The calculated density of the solution in (units)
    """
    # ensure ranges are okay here
    salt_to_id = {"None": 0, "hcl": 1, "nacl": 2, "naoh": 3}
    if salt.lower() in salt_to_id:
        return model.predict([[salt_to_id[salt.lower()], weight_percent*100, temperature]])[0]
    else:
        print("Pick a REAL salt, like a man.")