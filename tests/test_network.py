import json
from antiviral_qsar.network import Network


def test_snapshot_evidence_remains_immutable(config,httpx_mock):
    net=Network(config)
    url="https://example.org/source"
    httpx_mock.add_response(url=url,json={"value":1})
    first,path1=net.json(url,max_age=0)
    httpx_mock.add_response(url=url,json={"value":2})
    second,path2=net.json(url,max_age=0)
    assert path1 != path2
    assert json.loads((config.root/path1).read_bytes()) == {"value":1}
    cached,path3=net.json(url)
    assert cached == second and path3 == path2
    net.close()
