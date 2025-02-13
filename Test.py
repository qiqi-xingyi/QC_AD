# --*-- conding:utf-8 --*--
# @Time : 2/12/25 10:41 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Test.py

from qiskit_ibm_runtime import QiskitRuntimeService

if __name__ == '__main__':

    try:

        service = QiskitRuntimeService(
            channel='ibm_quantum',
            instance='ibm-q-ccf/qradle-catalyzer/qradle-catalyzer',
            token=''
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
