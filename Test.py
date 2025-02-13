# --*-- conding:utf-8 --*--
# @Time : 2/12/25 10:41 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Test.py

from qiskit_ibm_runtime import QiskitRuntimeService

def read_config(file_path):

    config = {}
    try:
        with open(file_path, "r") as file:
            for line in file:
                key, value = line.strip().split("=")
                config[key.strip()] = value.strip()
    except Exception as e:
        print(f"Fail to read token file: {e}")
        return None
    return config

if __name__ == '__main__':

    config_path = "config.txt"

    config = read_config(config_path)

    if not config or "TOKEN" not in config or "INSTANCE" not in config:
        print("The configuration file is missing a required field (TOKEN, INSTANCE)")
    else:

        try:

            service = QiskitRuntimeService(
                channel='ibm_quantum',
                instance=config["INSTANCE"],
                token=config["TOKEN"]
            )

            print("Connection successful: IBM Quantum account is active.")

            # Retrieve available quantum backends
            backends = service.backends(operational=True)
            if backends:
                print("Available IBM Quantum devices:")
                for backend in backends:
                    print(
                        f"  - {backend.name} (Status: {'Available' if backend.status().operational else 'Unavailable'})")
            else:
                print("No available quantum devices at the moment.")

        except Exception as e:
            print("Connection failed: Please check your API Key and network connection.")
            print("Error details:", e)
