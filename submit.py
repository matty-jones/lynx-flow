import subprocess
import time as T


def get_label_progress():
    status = subprocess.check_output(["python", "project.py", "status"]).decode("utf-8")
    # Split by line and trim empty lines
    status = [x for x in status.split("\n") if len(x) != 0]
    number_of_job_line = ["Total", "#", "of", "jobs:"]
    number_of_jobs = None
    recording = False
    label_data = {}
    for line in status:
        split_line = line.split()
        # Skip separator lines (that begin with a hyphen)
        if split_line[0][0] == "-":
            continue
        # Labels follow after 'label', 'progress'
        if ("label" in split_line) and ("progress" in split_line):
            recording = True
            continue
        if recording is True:
            # Everything beyond this point is a label
            label_name = split_line[0]
            label_percentage = float(split_line[-1][:-1])
            label_data[label_name] = label_percentage
    return label_data


if __name__ == "__main__":
    first_run = False
    second_run = False
    log_file_name = "./submit_output.log"
    # Reset the submit log
    with open(log_file_name, "w+") as log_file:
        pass
    while True:
        label_progress = get_label_progress()
        if (first_run is False) and ("generated" not in label_progress.keys()):
            print(
                "No simulations have been generated. Performing first run of project.py..."
            )
            first_run = True
        elif "generated" not in label_progress.keys():
            print(
                "`Generated' label not found. Probably nothing has finished yet. Sleeping..."
            )
            T.sleep(300)
        elif label_progress["generated"] != 100.0:
            print(
                "Current generation progress = {}%. Sleeping...".format(
                    label_progress["generated"]
                )
            )
            T.sleep(300)
        elif (
            (label_progress["generated"] == 100.0)
            and (second_run is False)
            and (label_progress["simulated"] != 100.0)
        ):
            print("Generate complete!")
            print(
                "All simulations have been generated, but not simulated. Performing second run of project.py..."
            )
            second_run = True
        elif (label_progress["generated"] == 100.0) and (
            label_progress["simulated"] != 100.0
        ):
            print(
                "Current simulation progress = {}%. Sleeping...".format(
                    label_progress["simulated"]
                )
            )
            T.sleep(300)
        elif (label_progress["generated"] == 100.0) and (
            label_progress["simulated"] == 100.0
        ):
            print("All generations and simulations completed! Exiting...")
            break
        with open(log_file_name, "a+") as log_file:
            subprocess.Popen(["python", "project.py", "submit"], stdout=log_file)
