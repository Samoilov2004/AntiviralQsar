import io
import json
import pytest
from rich.console import Console
from antiviral_qsar.progress import RunProgress, RemainingColumn


def test_resumed_progress_counts_saved_work_and_estimates_eta(tmp_path):
    clock = [0.0]
    terminal = Console(file=io.StringIO(), force_terminal=True, force_interactive=True, width=120)
    ui = RunProgress(terminal)
    ui.progress.get_time = lambda:clock[0]
    ui.log_path = tmp_path/'progress.log'
    with ui:
        for item in ui.track(range(5), 'Assays', completed=7):
            if item == 3:
                task = ui.progress.tasks[0]
                assert task.completed == 10 and task.total == 12
                assert RemainingColumn().render(task).plain == '2 left'
                assert task.time_remaining == 2
                snapshot = io.StringIO()
                Console(file=snapshot, width=120, color_system=None).print(ui.progress.get_renderable())
                rendered = snapshot.getvalue()
            clock[0] += 1
    assert 'ETA' in rendered and '10/12' in rendered and '2 left' in rendered
    entry = json.loads(ui.log_path.read_text().splitlines()[-1])
    assert entry['completed'] == 12 and entry['message'] == 'phase_finished'


def test_redirected_output_does_not_print_every_assay(tmp_path):
    output = io.StringIO()
    ui = RunProgress(Console(file=output, force_terminal=False, width=120))
    ui.log_path = tmp_path/'progress.log'
    for item in ui.track(range(1000), 'Assays'):
        ui.log('assay_started', id=f'CHEMBL{item}')
    assert len(output.getvalue().splitlines()) == 2
    assert 'CHEMBL' not in output.getvalue()
    assert '1,000/1,000 complete' in output.getvalue()
    assert 'CHEMBL999' in ui.log_path.read_text()


def test_interrupt_does_not_advance_unfinished_item(tmp_path):
    ui = RunProgress(Console(file=io.StringIO(), force_terminal=False))
    ui.log_path = tmp_path/'progress.log'
    with pytest.raises(KeyboardInterrupt):
        with ui:
            for item in ui.track(range(3), 'Assays'):
                if item == 1:
                    raise KeyboardInterrupt()
    entry = json.loads(ui.log_path.read_text().splitlines()[-1])
    assert entry['completed'] == 1 and entry['message'] == 'phase_interrupted'
    assert ui._depth == 0 and ui._stream is None


def test_paper_and_fragment_bars_share_one_live_display():
    ui = RunProgress(Console(file=io.StringIO(), force_terminal=True, width=100))
    with ui:
        for _ in ui.track(['paper'], 'Papers'):
            for _ in ui.track(['row'], 'Fragments'):
                assert [t.description for t in ui.progress.tasks] == ['Papers','Fragments']
    assert not ui.progress.tasks
