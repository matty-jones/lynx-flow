import flow
import environment


class LynxProject(flow.FlowProject):
    def __init__(self, *args, **kwargs):
        super(LynxProject, self).__init__(*args, **kwargs)
        env = flow.get_environment()
        self.add_operation(
            name='auto',
            cmd=lambda job: "python -u operations.py auto {}".format(job),
            post=[LynxProject.generated]
        )
        #self.add_operation(
        #    name='generate',
        #    cmd=lambda job: "python -u operations.py generate {}".format(job),
        #    post=[LynxProject.generated]
        #)
        #self.add_operation(
        #    name='simulate',
        #    cmd=lambda job: "python -u operations.py simulate {}".format(job),
        #    pre=[LynxProject.generated],
        #    post=[LynxProject.simulated]
        #)

    @flow.staticlabel()
    def generated(job):
        return job.isfile('output.hoomdxml')

    @flow.staticlabel()
    def simulated(job):
        return job.isfile('output_final.gsd')


if __name__ == '__main__':
    LynxProject().main()
