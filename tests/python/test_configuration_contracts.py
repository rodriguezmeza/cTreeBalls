"""Exercise public configuration entry points, precedence, diagnostics and recovery."""
from pathlib import Path
import re
import subprocess
import pytest
from cyballs import cballs

ROOT = Path(__file__).resolve().parents[2]
INVALID = ['truejunk', 'falsejunk', 'garbage', '', 'yes please', '2', '-1']
VALID = [('TRUE', True), ('t', True), ('YeS', True), ('Y', True), ('1', True),
         ('false', False), ('F', False), ('NO', False), ('n', False), ('0', False)]


def params(path):
    return dict(searchMethod='octree-2balls-omp', nbody=8, numberThreads=1,
                rootDir=str(path), verbose=0, verbose_log=0, usePeriodic=False,
                useLogHist=False, rangeN=1, rminHist=.01, sizeHistN=4,
                mChebyshev=2, options='no-out-Hist,no-smooth-pivot')


def cli(parameters, *args):
    return subprocess.run([str(ROOT/'cballs'), *args,
                           *[f'{k}={v}' for k,v in parameters.items()]],
                          text=True, capture_output=True, timeout=30)


@pytest.mark.parametrize('key', ['usePeriodic', 'useLogHist'])
@pytest.mark.parametrize('token', INVALID)
def test_cli_rejects_entire_malformed_boolean(tmp_path, key, token):
    p=params(tmp_path);p[key]=token
    result=cli(p)
    assert result.returncode != 0
    assert f'invalid boolean value for {key}' in result.stdout+result.stderr


@pytest.mark.parametrize('key', ['usePeriodic', 'useLogHist'])
@pytest.mark.parametrize('token', INVALID)
def test_python_rejects_and_same_object_recovers(tmp_path, key, token):
    model=cballs()
    try:
        p=params(tmp_path);p[key]=token;model.set(p)
        with pytest.raises(Exception, match=f'invalid boolean value for {key}'):
            model.Run()
        assert not model.state and not model.getBodytableAllocated()
        model.set(**{key: False});model.Run()
        assert model.state
    finally:
        model.struct_cleanup()


@pytest.mark.parametrize('key', ['usePeriodic', 'useLogHist'])
@pytest.mark.parametrize('token,expected', VALID)
def test_complete_boolean_aliases_preserve_value(tmp_path, key, token, expected):
    model=cballs()
    try:
        p=params(tmp_path);p[key]=token;model.set(p);model.Run()
        assert model.run_settings['effective'][key] is expected
    finally:
        model.struct_cleanup()


def parameter_file(tmp_path):
    result=cli(params(tmp_path/'seed'))
    assert result.returncode == 0, result.stdout+result.stderr
    return (tmp_path/'seed/parameters_null-cballs-usedvalues').read_text()


@pytest.mark.parametrize('key,token', [('usePeriodic','truejunk'), ('useLogHist','falsejunk')])
def test_file_rejects_boolean_with_key(tmp_path, key, token):
    text=parameter_file(tmp_path)
    text=re.sub(r'^'+key+r'\s*=.*$', key+' = '+token, text, flags=re.M)
    path=tmp_path/'bad.params';path.write_text(text)
    result=cli({},str(path))
    assert result.returncode != 0
    assert f'invalid boolean value for {key}' in result.stdout+result.stderr


@pytest.mark.parametrize('declaration,canonical', [('sizeHistN = 7', 'sizeHistN'),
                                                   ('search = octree-2balls-omp', 'searchMethod')])
def test_duplicate_file_keys_including_aliases_have_line_numbers(tmp_path, declaration, canonical):
    text=parameter_file(tmp_path);path=tmp_path/'duplicate.params'
    path.write_text(text+'\n'+declaration+'\n')
    result=cli({},str(path))
    message=result.stdout+result.stderr
    assert result.returncode != 0
    assert 'unknown or duplicated' in message and canonical in message
    assert re.search(r'duplicate\.params:\d+:',message)
    assert re.search(r'first declared at line \d+',message)


def test_cli_override_is_allowed_after_unique_file(tmp_path):
    path=tmp_path/'valid.params';path.write_text(parameter_file(tmp_path))
    result=cli({'sizeHistN':7,'rootDir':tmp_path/'override'},str(path))
    assert result.returncode == 0, result.stdout+result.stderr
    used=(tmp_path/'override/parameters_null-cballs-usedvalues').read_text()
    assert re.search(r'^sizeHistN\s*=\s*7\s*$',used,re.M)


def native_mpi_checks(directory, launcher):
    """A bad rank must terminate collectively; corrected files then run successfully."""
    import json
    import shlex
    import sys
    directory.mkdir(parents=True, exist_ok=False)
    template=parameter_file(directory)
    good=directory/'valid.params';good.write_text(template)
    variants={
        'boolean': re.sub(r'^usePeriodic\s*=.*$', 'usePeriodic = truejunk', template, flags=re.M),
        'duplicate': template+'\nsizeHistN = 7\n',
        'alias': template+'\nsearch = octree-2balls-omp\n',
    }
    records=[]
    for label,text in variants.items():
        bad=directory/(label+'.params');bad.write_text(text)
        for rank in (0,1):
            for failure in (True,False):
                argv=shlex.split(launcher)+['-n','2',sys.executable,str(Path(__file__).resolve()),
                      '--native-mpi-worker',str(bad if failure else good),str(good),str(rank)]
                result=subprocess.run(argv,text=True,capture_output=True,timeout=60)
                output=result.stdout+result.stderr
                log=f'{label}-rank{rank}-'+('failure' if failure else 'recovery')+'.log'
                (directory/log).write_text(output)
                if failure:
                    expected='invalid boolean value' if label=='boolean' else 'unknown or duplicated'
                    if result.returncode==0 or expected not in output:
                        raise RuntimeError(f'MPI configuration failure was not propagated: {log}')
                elif result.returncode:
                    raise RuntimeError(f'MPI configuration recovery failed: {log}')
                records.append(dict(case=label,rank=rank,failure=failure,returncode=result.returncode,log=log))
    (directory/'results.json').write_text(json.dumps(dict(status='PASS',cases=records),indent=2)+'\n')
    print('PASS: 12 native MPI configuration failure/retry executions')


if __name__ == '__main__':
    import argparse
    import os
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='--native-mpi-worker':
        rank=int(next(os.environ[k] for k in ('OMPI_COMM_WORLD_RANK','PMI_RANK','PMIX_RANK') if k in os.environ))
        path=sys.argv[2] if rank==int(sys.argv[4]) else sys.argv[3]
        binary=str(ROOT/'cballs')
        os.execv(binary,[binary,path,'searchMethod=octree-2balls-mpi','seed=123'])
    parser=argparse.ArgumentParser()
    parser.add_argument('--mpi-command',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    native_mpi_checks(args.output,args.mpi_command)


@pytest.mark.parametrize('key', ['verbose', 'verbose_log'])
@pytest.mark.parametrize('value', [32768, -32769])
def test_short_parameter_range_and_recovery(tmp_path, key, value):
    p=params(tmp_path);p[key]=value
    result=cli(p)
    assert result.returncode != 0 and 'outside the short integer range' in result.stdout+result.stderr
    model=cballs()
    try:
        model.set(p)
        with pytest.raises(Exception, match='outside the short integer range'):
            model.Run()
        model.set(**{key: 0});model.Run()
        assert model.state
    finally:
        model.struct_cleanup()
