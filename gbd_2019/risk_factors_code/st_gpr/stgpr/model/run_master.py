import argparse
import dataclasses
import sys

import orm_stgpr
from orm_stgpr.lib.constants import parameters
from orm_stgpr.lib.util import old

from stgpr.model import jobmon_workflow


@dataclasses.dataclass
class MasterArgs:
    run_id: int
    error_log_path: str
    output_log_path: str
    project: str
    nparallel: int
    code_version: str
    output_path: str


def parse_args() -> MasterArgs:
    parser = argparse.ArgumentParser()
    parser.add_argument('run_id', type=int)
    parser.add_argument('error_log_path', type=str)
    parser.add_argument('output_log_path', type=str)
    parser.add_argument('project', type=str)
    parser.add_argument('nparallel', type=int)
    parser.add_argument('code_version', type=str)
    parser.add_argument('output_path', type=str)
    args = vars(parser.parse_args())
    return MasterArgs(**args)


def run_job_swarm(
        run_id: int,
        error_log_path: str,
        output_log_path: str,
        project: str,
        nparallel: int,
        code_version: str,
        output_path: str
) -> None:
    params = orm_stgpr.get_parameters(run_id)
    n_param_sets = old.determine_n_parameter_sets(params)
    swarm = jobmon_workflow.STGPRJobSwarm(
        run_id=run_id,
        run_type=params[parameters.RUN_TYPE],
        holdouts=params[parameters.HOLDOUTS],
        draws=params[parameters.GPR_DRAWS],
        nparallel=nparallel,
        n_parameter_sets=n_param_sets,
        cluster_project=project,
        error_log_path=error_log_path,
        output_log_path=output_log_path,
        location_set_id=params[parameters.LOCATION_SET_ID],
        gbd_round_id=params[parameters.GBD_ROUND_ID],
        custom_stage1=params[parameters.CUSTOM_STAGE_1],
        rake_logit=params[parameters.RAKE_LOGIT],
        code_version=code_version,
        decomp_step=params[parameters.DECOMP_STEP],
        modelable_entity_id=params[parameters.MODELABLE_ENTITY_ID],
        output_path=output_path
    )
    status = swarm.run()
    sys.exit(status)


def main() -> None:
    args = parse_args()
    run_job_swarm(
        args.run_id,
        args.error_log_path,
        args.output_log_path,
        args.project,
        args.nparallel,
        args.code_version,
        args.output_path
    )


if __name__ == '__main__':
    main()
